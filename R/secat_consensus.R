# ==============================================================================
# SCRIPT:   secat_consensus.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# PURPOSE:  Compute the maximal overlapping 16S consensus region across studies
#
# OVERVIEW:
#   Identifies the "Consensus Region" -- the shared genomic window that allows
#   disparate studies targeting different 16S variable regions to be compared.
#   Two strategies are implemented: (1) Primer Mode aligns physical primer
#   sequences to SILVA to determine canonical amplicon coordinates, then takes
#   the intersection; (2) Study Mode builds a compatibility graph from empirical
#   study coordinates and finds the largest clique of mutually overlapping
#   studies, with greedy optimisation to remove outlier "constrictor" studies.
#
# SECTIONS:
#   1. Primer Mode (Physical Alignment) -- vmatchPattern-based coordinate calling
#   2. Study Mode (Graph-Theoretic Clique Detection) -- adjacency matrix + greedy
#      constrictor removal via recursive optimisation
#
# SOURCED BY:
#   - Nextflow processes that call consensus region calculation
#   - 05_calculate_consensus.R (or equivalent orchestration script)
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: calculate_mode
# Purpose:  Statistical mode (most frequent value) of a vector
#
# Parameters:
#   @param x [vector] - Numeric or character vector
#
# Returns:
#   @return [any] - Most frequent value; NA if empty. Ties broken by first max.
#
# Notes:
#   Duplicated from secat_utils.R to keep this file self-contained.
# ------------------------------------------------------------------------------
calculate_mode <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ------------------------------------------------------------------------------
# Function: relax_consensus_coords
# Purpose:  Validate and fix inverted consensus coordinates (start >= end)
#
# Parameters:
#   @param clique_result [list] - Output from find_largest_overlapping_clique()
#   @param all_starts    [numeric] - Start coordinates for all studies
#   @param all_ends      [numeric] - End coordinates for all studies
#   @param study_names   [character] - Study identifiers
#   @param study_name    [character] - Label for log messages
#
# Returns:
#   @return [list] - Corrected clique_result with guaranteed start < end
#
# Notes:
#   Fallback uses the union of all study bounds with a minimum 50 bp window.
#   This prevents downstream crashes from invalid coordinate ranges.
# ------------------------------------------------------------------------------
relax_consensus_coords <- function(clique_result, all_starts, all_ends, study_names, study_name = "consensus") {
  cons_start <- clique_result$start
  cons_end <- clique_result$end

  if (!is.na(cons_start) && !is.na(cons_end) && cons_start >= cons_end) {
    message(sprintf("[WARN] Invalid Consensus Detected (Start: %d >= End: %d) for %s. Applying smart relaxation...",
                    cons_start, cons_end, study_name))

    # Fallback: use the full union of all study bounds
    study_min <- min(all_starts, na.rm = TRUE)
    study_max <- max(all_ends, na.rm = TRUE)

    # Guarantee at least 50 bp window (matches min_overlap default)
    relaxed_start <- min(study_min, study_max - 50)
    relaxed_end <- max(study_max, study_min + 50)

    clique_result$start <- relaxed_start
    clique_result$end <- relaxed_end

    message(sprintf("  -> Relaxed to study bounds [%d,%d] (length: %d bp)",
                    relaxed_start, relaxed_end, relaxed_end - relaxed_start))
  }
  return(clique_result)
}

# ==== SECTION 1: PRIMER MODE (Physical Alignment) ====
# Aligns physical primer sequences to SILVA reference to determine canonical
# amplicon coordinates. The consensus is the intersection of all amplicons
# (max start, min end).

# ------------------------------------------------------------------------------
# Function: find_consensus_region
# Purpose:  Determine consensus region by aligning primer sequences to reference
#
# Parameters:
#   @param primer_info_df   [data.frame] - Columns: primer_name, fwd_seq, rev_seq
#   @param reference_db_path [character] - Path to reference FASTA (SILVA RNA/DNA)
#   @param config            [list]      - Optional config overrides
#
# Returns:
#   @return [data.frame] - Primer coordinates with global consensus_start/end columns
#
# Notes:
#   Uses Biostrings::vmatchPattern with 4 mismatches to locate each primer in
#   the ungapped reference. Modal hit positions are used (robust to paralogs).
#   The global consensus = max(all starts) to min(all ends), i.e. the region
#   covered by every primer pair simultaneously. A pre-ungapped reference can
#   be passed via SECAT_UNGAPPED_DB env var to skip in-memory gap removal.
# ------------------------------------------------------------------------------
find_consensus_region <- function(primer_info_df, reference_db_path, config = list()) {

  # --- 1. Input Validation ---
  if (!file.exists(reference_db_path)) {
    stop("Error: Reference database file not found at: ", reference_db_path)
  }

 # --- 2. Load Reference Database ---
  # Check if a pre-computed ungapped DB path was provided by the shell (Phase 1 optimization)
  ungapped_db_path <- Sys.getenv("SECAT_UNGAPPED_DB")
  
  if (ungapped_db_path != "" && file.exists(ungapped_db_path)) {
      # FAST PATH: Load pre-cleaned DB from disk (created by sed in SGE script)
      message("Loading pre-ungapped reference from temp file (fast)...")
      reference_sequences_ungapped <- Biostrings::readDNAStringSet(ungapped_db_path)
      message(paste("...database loaded. Working with", length(reference_sequences_ungapped), "sequences."))
      
  } else {
      # FALLBACK PATH: In-memory processing (Memory intensive, used if SGE script didn't run sed)
      message("Loading reference database (as RNA)...")
      reference_sequences_raw <- Biostrings::readRNAStringSet(reference_db_path)
      
      message("Converting RNA database to DNA...")
      reference_sequences <- Biostrings::DNAStringSet(reference_sequences_raw)
      
      # Clean up raw RNA object
      rm(reference_sequences_raw)
      gc()
      
      message("Removing gaps (in-memory fallback)...")
      # Use robust character-based replacement instead of function calls that might not export
      reference_sequences_ungapped <- Biostrings::DNAStringSet(gsub("-", "", as.character(reference_sequences)))
      
      # Clean up gapped DNA object
      rm(reference_sequences)
      gc()
      
      message(paste("...gaps removed. Working with", length(reference_sequences_ungapped), "sequences."))
  }
  
  # --- 4. Initialize results table ---
  primer_coordinates <- tibble::tibble(
    primer_name = character(),
    primer_start = integer(),
    primer_end = integer()
  )

  # --- 5. Main Loop: Process each unique primer set ---
  unique_primers <- primer_info_df %>% dplyr::distinct(primer_name, fwd_seq, rev_seq)

  message("Finding primer coordinates for each unique primer set...")
  for (i in 1:nrow(unique_primers)) {
    primer_name <- unique_primers$primer_name[i]
    message(paste0("  Aligning ", primer_name, "..."))

    # Sanitize primer sequences (remove non-IUPAC characters)
    fwd_seq_sanitized <- gsub("[^ACGTMRWSYKVHDBN]", "", toupper(unique_primers$fwd_seq[i]))
    rev_seq_sanitized <- gsub("[^ACGTMRWSYKVHDBN]", "", toupper(unique_primers$rev_seq[i]))

    # Use ungapped reference for vmatchPattern
    # Allowing 4 mismatches accounts for degenerate bases and minor biological variation.
    fwd_hits <- Biostrings::vmatchPattern(
      Biostrings::DNAString(fwd_seq_sanitized),
      reference_sequences_ungapped,
      max.mismatch = 4
    )
    rev_hits <- Biostrings::vmatchPattern(
      Biostrings::reverseComplement(Biostrings::DNAString(rev_seq_sanitized)),
      reference_sequences_ungapped,
      max.mismatch = 4
    )

    # Use the mode (most common) position to define the canonical coordinate.
    modal_start <- calculate_mode(unlist(Biostrings::startIndex(fwd_hits)))
    modal_end <- calculate_mode(unlist(Biostrings::endIndex(rev_hits)))

    if (is.na(modal_start) || is.na(modal_end)) {
      warning(paste("Warning: Primers for '", primer_name, "' not found in the reference database. Skipping."))
      next
    }

    primer_coordinates <- primer_coordinates %>%
      tibble::add_row(primer_name = primer_name, primer_start = modal_start, primer_end = modal_end)
  }

  # --- 6. Calculate global consensus region ---
  # The intersection of ALL valid amplicons.
  # Start = The right-most (maximum) start position.
  # End   = The left-most (minimum) end position.
  message("Calculating global consensus region...")
  if (nrow(primer_coordinates) == 0) {
    warning("No primers were successfully mapped. Cannot calculate consensus region.")
    return(primer_coordinates)
  }

  consensus_start_val <- max(primer_coordinates$primer_start)
  consensus_end_val   <- min(primer_coordinates$primer_end)
  message(paste("Global consensus region determined to be:", consensus_start_val, "-", consensus_end_val))

  primer_coordinates_final <- primer_coordinates %>%
    dplyr::mutate(
      consensus_start = consensus_start_val,
      consensus_end = consensus_end_val
    )

  return(primer_coordinates_final)
}

#===============================================================================
# SECTION 2: STUDY MODE (Graph-Theoretic Clique Detection)
#===============================================================================

#-------------------------------------------------------------------------------
# Function:   find_largest_overlapping_clique
#-------------------------------------------------------------------------------
# Description:
#   Identifies the largest set of studies that can form a valid meta-analysis
#   group (a "clique") based on coordinate overlap.
#
# Scientific Context:
#   In "Study Mode", we have many datasets with potentially disparate coordinates.
#   We want to include as many studies as possible in the meta-analysis, but
#   they MUST all overlap by at least `min_overlap` bp. This creates a trade-off:
#   including an outlier study might force the consensus window to be tiny.
#   This function finds the optimal balance.
#
# Algorithm:
#   1. Adjacency Matrix: Build a graph where nodes are studies and edges exist
#      if the pair overlaps by >= `min_overlap`.
#   2. Initial Clique: Find the largest subset of mutually overlapping studies.
#   3. Optimization: Calls `optimize_clique_recursive` to refine this set. It
#      checks if dropping "constrictor" studies (those at the edges) yields a
#      significant gain in consensus length.
#
# Parameters:
#   @param starts [vector] - Start coordinates for each study.
#   @param ends [vector] - End coordinates for each study.
#   @param study_names [vector] - IDs/names of the studies.
#   @param min_overlap [integer] - Minimum bp required to consider an overlap valid.
#
# Returns:
#   @return [list] - Metadata including the final start/end, list of included
#           studies, and list of excluded studies.
#-------------------------------------------------------------------------------
find_largest_overlapping_clique <- function(starts, ends, study_names, min_overlap = 50) {

  n <- length(starts)
  if (n == 0) return(NULL)

  # --- PHASE 1: Build Compatibility Graph ---
  # A pair is compatible if their intersection is at least `min_overlap` long.
  adjacency <- matrix(FALSE, nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      inter_start <- max(starts[i], starts[j])
      inter_end <- min(ends[i], ends[j])
      if ((inter_end - inter_start) >= min_overlap) {
        adjacency[i, j] <- TRUE
      }
    }
  }

  # --- PHASE 2: Find Largest Initial Clique ---
  # Simple greedy search for the maximum clique.
  # (Note: Max Clique is NP-hard, but for N < 1000 studies, this heuristic is fast enough).
  best_clique <- c()

  for (i in 1:n) {
    candidate <- c(i)
    for (j in 1:n) {
      if (i == j) next
      # If study J is connected to ALL studies currently in the candidate set...
      if (all(adjacency[j, candidate])) {
        candidate <- c(candidate, j)
      }
    }
    if (length(candidate) > length(best_clique)) best_clique <- candidate
  }

  # --- PHASE 3: Iterative Optimization (Remove Constrictors) ---
  # Load threshold from config (if available) or default to 20% gain.
  # This means: "Only drop a study if it increases consensus length by >20%".
  threshold <- if(exists("CONSENSUS_OPTIMIZATION_THRESHOLD")) CONSENSUS_OPTIMIZATION_THRESHOLD else 0.20
  min_studies <- if(exists("MIN_CONSENSUS_STUDIES")) MIN_CONSENSUS_STUDIES else 3

  optimized_clique <- optimize_clique_recursive(
    indices = best_clique,
    starts = starts,
    ends = ends,
    threshold = threshold,
    min_studies = min_studies
  )

# --- PHASE 4: Final Output with Validation ---
final_indices <- optimized_clique
if (length(final_indices) == 0) final_indices <- best_clique

clique_result <- list(
  start = max(starts[final_indices]),
  end = min(ends[final_indices]),
  n_studies = length(final_indices),
  included_studies = study_names[final_indices],
  excluded_studies = study_names[-final_indices]
)

# CRITICAL: Apply relaxation BEFORE returning to CSV
clique_result <- relax_consensus_coords(clique_result, starts, ends, study_names, "Family")
return(clique_result)
}

#-------------------------------------------------------------------------------
# Function:   optimize_clique_recursive
#-------------------------------------------------------------------------------
# Description:
#   Recursively removes "constrictor" studies to optimize the consensus window.
#
# Algorithm:
#   1. Identify current consensus length.
#   2. Identify the "Start Constrictor" (study with the highest start position).
#   3. Identify the "End Constrictor" (study with the lowest end position).
#   4. Simulate removing each. Calculate the % gain in consensus length.
#   5. If the best gain exceeds `threshold`, commit the removal and recurse.
#
# Parameters:
#   @param indices [vector] - Indices of studies currently in the clique.
#   @param threshold [numeric] - % improvement required to drop a study (0.20 = 20%).
#
# Returns:
#   @return [vector] - The optimized list of indices.
#-------------------------------------------------------------------------------
optimize_clique_recursive <- function(indices, starts, ends, threshold, min_studies) {

  # Stop recursion if we hit minimum study count (safety floor).
  if (length(indices) <= min_studies) return(indices)

  # 1. Calculate Current State
  curr_starts <- starts[indices]
  curr_ends <- ends[indices]

  current_consensus_start <- max(curr_starts)
  current_consensus_end   <- min(curr_ends)
  current_len <- current_consensus_end - current_consensus_start

  # If current length is negative/zero (invalid), any valid length is infinite improvement
  if (current_len <= 0) current_len <- 1

  # 2. Identify Constrictors (The studies defining the boundaries)
  bad_start_local_idx <- which.max(curr_starts)
  bad_end_local_idx <- which.min(curr_ends)

  # 3. Simulate Removing "Start Constrictor"
  test_indices_start <- indices[-bad_start_local_idx]
  len_drop_start <- (min(ends[test_indices_start]) - max(starts[test_indices_start]))
  gain_start <- (len_drop_start - current_len) / current_len

  # 4. Simulate Removing "End Constrictor"
  test_indices_end <- indices[-bad_end_local_idx]
  len_drop_end <- (min(ends[test_indices_end]) - max(starts[test_indices_end]))
  gain_end <- (len_drop_end - current_len) / current_len

  # 5. Evaluate and Recurse
  best_gain <- max(gain_start, gain_end)

  if (best_gain > threshold) {
    # Found significant improvement
    if (gain_start >= gain_end) {
      return(optimize_clique_recursive(test_indices_start, starts, ends, threshold, min_studies))
    } else {
      return(optimize_clique_recursive(test_indices_end, starts, ends, threshold, min_studies))
    }
  }

  # No significant improvement, return current set
  return(indices)
}
