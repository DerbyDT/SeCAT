# ==============================================================================
# SCRIPT:   secat_mapping.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# PURPOSE:  Reference mapping engine — aligns study ASV sequences to the SILVA
#           reference to determine each study's amplicon coordinates.
#
# OVERVIEW:
#   This library provides the core function map_study_to_reference(), which:
#   1. Optionally subsamples ASVs for speed (configurable via USE_ALL_ASVS)
#   2. Builds or loads a SILVA reference subset for DECIPHER alignment
#   3. Aligns study ASVs to the SILVA profile using DECIPHER::AlignSeqs()
#   4. Determines the modal start/end positions across all aligned ASVs
#   5. Returns the study-level coordinates, per-ASV coordinates, and alignment
#
#   The modal coordinate approach is robust to outlier ASVs that align poorly
#   (e.g., chimeric sequences or non-16S contaminants). The alignment is also
#   saved for use in the trimming analysis (06_analyse_real.R), where sequences
#   are trimmed in alignment space rather than raw sequence space.
#
# SOURCED BY:
#   - 02_study_mapping.R (via map_study_to_reference())
# ==============================================================================

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(tibble))

# ------------------------------------------------------------------------------
# Function: calculate_mode
# ------------------------------------------------------------------------------
# Purpose:
#   Returns the statistical mode (most frequent value) of a numeric vector.
#   Used to find the dominant start/end alignment column across all ASVs in a
#   study, which defines the consensus amplicon window.
#
# Parameters:
#   @param x  [numeric] Vector of values (typically alignment column indices).
#
# Returns:
#   @return [numeric(1)] The most frequent value, or NA_real_ if input is empty.
#
# Notes:
#   Ties are broken by first occurrence (which.max returns the first maximum).
#   This is acceptable because tied modes in amplicon data are rare — the vast
#   majority of ASVs from a single primer pair share identical start/end columns.
# ------------------------------------------------------------------------------
calculate_mode <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ------------------------------------------------------------------------------
# Function: parse_primer_positions (Legacy Helper)
# ------------------------------------------------------------------------------
# Purpose:
#   Extracts numeric start/end coordinates from standardised primer name strings
#   (e.g., "515F_806R" -> start=515, end=806). This is the legacy/fallback
#   approach that trusts primer naming conventions rather than empirical alignment.
#
# Parameters:
#   @param primer_name [character(1)] Primer pair string with forward (F) and
#     reverse (R) tokens separated by "_", "-", or ".".
#
# Returns:
#   @return [list] With elements $start and $end (integer or NA_integer_).
#
# Notes:
#   Only used in "primer" analysis mode. The "study" mode (DECIPHER alignment)
#   is preferred because primer names can be non-standard or misleading.
# ------------------------------------------------------------------------------
parse_primer_positions <- function(primer_name) {
    # Allow separators like "_", "-", or "." and extra prefixes
    parts <- strsplit(primer_name, "[-_.]")[[1]]
    parts <- parts[nzchar(parts)]

    # Find the first forward- and reverse-like tokens
    fwd_token <- parts[grep("F$", parts, ignore.case = TRUE)[1]]
    rev_token <- parts[grep("R$", parts, ignore.case = TRUE)[1]]

    get_pos <- function(x) {
        if (is.na(x)) return(NA_integer_)
        m <- regmatches(x, regexpr("[0-9]+", x))
        if (length(m) > 0) as.integer(m) else NA_integer_
    }

    fwd_pos <- get_pos(fwd_token)
    rev_pos <- get_pos(rev_token)

    list(start = fwd_pos, end = rev_pos)
}

# ------------------------------------------------------------------------------
# Function: calculate_consensus_coordinates
# ------------------------------------------------------------------------------
# Purpose:
#   Aggregates per-ASV alignment coordinates into a single consensus amplicon
#   window for the study. This is the critical step that converts hundreds or
#   thousands of individual ASV mappings into a single (start, end) pair.
#
# Parameters:
#   @param all_ref_starts [integer] Vector of alignment-column start positions
#     for each successfully mapped ASV.
#   @param all_ref_ends   [integer] Vector of alignment-column end positions.
#   @param method         [character(1)] Consensus strategy; one of:
#     - "modal"      (Default, recommended): Most frequent start/end values.
#                    Best for amplicon data where most ASVs share exact primer
#                    binding sites. Robust to chimeric or off-target outlier ASVs.
#     - "percentile": Uses 10th percentile for start, 90th for end. More
#                    conservative; captures slight primer-landing variation.
#     - "minimum":   min(start) to max(end). Maximises coverage but vulnerable
#                    to a single misaligned ASV widening the window.
#
# Returns:
#   @return [list] With elements $start and $end (integer).
#
# Notes:
#   The method is set via config$STUDY_ALIGNMENT_METHOD (default "modal").
#   For well-behaved amplicon datasets, all three methods yield identical results
#   because ASVs from a single primer pair align to the same reference columns.
# ------------------------------------------------------------------------------
calculate_consensus_coordinates <- function(all_ref_starts, all_ref_ends, method = "modal") {
  if (method == "modal") {
    start <- calculate_mode(all_ref_starts)
    end <- calculate_mode(all_ref_ends)
  } else if (method == "percentile") {
    start <- quantile(all_ref_starts, 0.10, type = 1)
    end <- quantile(all_ref_ends, 0.90, type = 1)
  } else if (method == "minimum") {
    start <- min(all_ref_starts)
    end <- max(all_ref_ends)
  } else {
    warning(paste("Unknown method:", method, ". Defaulting to 'modal'"))
    start <- calculate_mode(all_ref_starts)
    end <- calculate_mode(all_ref_ends)
  }
  return(list(start = start, end = end))
}

# ------------------------------------------------------------------------------
# Function: map_study_to_reference
# ------------------------------------------------------------------------------
# Purpose:
#   Core mapping function for the SeCAT pipeline. Determines the amplicon
#   coordinates of a study by aligning its ASV sequences to a SILVA reference
#   profile. This is the empirical alternative to trusting primer-name metadata.
#
# Scientific Context:
#   Accurate coordinate mapping is essential for cross-study meta-analysis. A
#   study labelled "V4" may actually cover V3-V4 due to a protocol modification;
#   relying on metadata alone would lead to invalid biological comparisons.
#   This function verifies the *actual* amplicon coverage empirically by aligning
#   real ASV sequences against a common reference coordinate system (SILVA).
#
# Algorithm (Study Mode — the recommended path):
#   1. Subsampling: Randomly selects N ASVs (default controlled by ASV_SAMPLE_SIZE
#      or USE_ALL_ASVS) to represent the study. A subset is statistically
#      sufficient to find the consensus window and dramatically reduces runtime.
#   2. Reference Prep: Loads the pre-aligned SILVA reference, optionally
#      subsetting it (REFERENCE_ALIGNMENT_MODE / REFERENCE_SUBSET_SIZE) to
#      accelerate profile construction.
#   3. Alignment: First aligns ASVs to each other via DECIPHER::AlignSeqs(),
#      then merges the ASV profile into the SILVA profile via
#      DECIPHER::AlignProfiles(). This preserves the reference coordinate system.
#   4. Coordinate Extraction: Finds the first and last non-gap character for
#      each aligned ASV, yielding per-ASV start/end column indices.
#   5. Consensus: Applies calculate_consensus_coordinates() (modal by default)
#      to collapse per-ASV coordinates into a single study-level window.
#
# Parameters:
#   @param study_info        [list] Study metadata; must contain $study_name,
#     $asv_fasta_path, $primer_name, $initial_fwd_trim, $initial_rev_trim.
#   @param reference_db_path [character(1)] Path to the SILVA aligned FASTA file.
#   @param config            [list] Pipeline configuration; key fields used:
#     $ANALYSIS_MODE             — "study" (DECIPHER) or "primer" (legacy).
#     $USE_ALL_ASVS              — logical; if TRUE, skip subsampling.
#     $ASV_SAMPLE_SIZE           — integer; max ASVs to subsample.
#     $REFERENCE_ALIGNMENT_MODE  — "full" or "subset".
#     $REFERENCE_SUBSET_SIZE     — integer; how many SILVA sequences to use.
#     $STUDY_ALIGNMENT_METHOD    — "modal", "percentile", or "minimum".
#
# Returns:
#   @return [list] with three elements:
#     $summary   — tibble: one row with study-level consensus coordinates,
#                  amplicon length, method used, and number of ASVs mapped.
#     $coords    — tibble: per-ASV start/end coordinates for diagnostics.
#     $alignment — DNAStringSet: the aligned ASV sequences in SILVA coordinate
#                  space. Used downstream in 06_analyse_real.R for trimming
#                  in alignment space.
#   Returns NULL if the FASTA is missing/empty or alignment fails.
# ------------------------------------------------------------------------------
map_study_to_reference <- function(study_info, reference_db_path, config) {

  message(paste("Processing study:", study_info$study_name))

  # Determine Mode from Config (Default to "study" if missing)
  mode <- if(!is.null(config$ANALYSIS_MODE)) config$ANALYSIS_MODE else "study"
  message(paste("  - Analysis Mode:", mode))

  # --- LEGACY PATH: PRIMER MODE ---
  # Trusts the primer name string to define coordinates.
  if (mode == "primer") {
      message("  - Using theoretical primer coordinates (Legacy Mode).")

      # 1. Parse primer positions
      primer_positions <- parse_primer_positions(study_info$primer_name)
      if (is.na(primer_positions$start) || is.na(primer_positions$end)) {
        warning(paste("Could not parse positions from primer name:", study_info$primer_name))
        return(NULL)
      }

      # 2. Return theoretical coordinates wrapped in list structure
      # We assume these apply directly to the reference genome/alignment
      summary_df <- tibble::tibble(
        study_name = study_info$study_name,
        primer_name = study_info$primer_name,
        ref_start_method = "PRIMER_THEORETICAL",
        ref_start = primer_positions$start,
        ref_end = primer_positions$end,
        analysis_amplicon_length = primer_positions$end - primer_positions$start + 1,
        original_start_theoretical = primer_positions$start,
        initial_fwd_trim = study_info$initial_fwd_trim,
        initial_rev_trim = study_info$initial_rev_trim,
        num_asvs_mapped = NA # Not mapping ASVs
      )

      return(list(
        summary = summary_df,
        coords = NULL,
        alignment = NULL
      ))
  }

  # --- MODERN PATH: STUDY/ALIGNMENT MODE (DECIPHER) ---
  # Empirically aligns ASVs to the reference profile.

  # 1. Check Dependencies
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    stop("Package 'DECIPHER' is required but not installed.")
  }

  # 2. Validate Study ASVs
  if (!file.exists(study_info$asv_fasta_path)) {
    warning(paste("FASTA not found for", study_info$study_name, ". Skipping."))
    return(NULL)
  }

  asv_sequences <- Biostrings::readDNAStringSet(study_info$asv_fasta_path)
  if (length(asv_sequences) == 0) {
    warning(paste("FASTA is empty for", study_info$study_name, ". Skipping."))
    return(NULL)
  }
  message(paste("  - Loaded", length(asv_sequences), "ASV sequences."))

  # 3. ASV Subsampling Logic
  # Aligning millions of ASVs is slow; a random subset defines the window accurately.
  if (config$USE_ALL_ASVS) {
     message("  - Using ALL ASVs for alignment (High Compute Mode).")
  } else {
     sample_size <- min(config$ASV_SAMPLE_SIZE, length(asv_sequences))
     if (length(asv_sequences) > sample_size) {
       set.seed(42)
       asv_sample_indices <- sample(seq_along(asv_sequences), sample_size)
       asv_sequences <- asv_sequences[asv_sample_indices]
       message(paste("  - Subsampled to", length(asv_sequences), "ASVs for mapping."))
     }
  }

  # 4. Robust Reference Loading
  message("  - Loading Reference Database...")
  # Using readBStringSet allows flexibility (works for RNA or DNA files)
  ref_raw_bstring <- Biostrings::readBStringSet(reference_db_path)

  # Optionally subset the reference to speed up Profile creation
  if (config$REFERENCE_ALIGNMENT_MODE == "full") {
      message("  - Using FULL reference database as profile.")
      ref_subset_bstring <- ref_raw_bstring
  } else {
      subset_size <- config$REFERENCE_SUBSET_SIZE
      actual_size <- min(subset_size, length(ref_raw_bstring))
      message(paste("  - Using SUBSET of reference database (", actual_size, "sequences) as profile."))
      set.seed(123)
      ref_indices <- sample(seq_along(ref_raw_bstring), actual_size)
      ref_subset_bstring <- ref_raw_bstring[ref_indices]
  }

  message("  - Standardizing reference profile format...")
  # Ensure standard DNA alphabet (T not U) and gap characters (-)
  clean_seqs <- function(bstring) {
      seqs_char <- as.character(bstring)
      seqs_char <- gsub(".", "-", seqs_char, fixed = TRUE)
      seqs_char <- gsub("U", "T", seqs_char, ignore.case = TRUE)
      return(seqs_char)
  }
  cleaned_chars <- clean_seqs(ref_subset_bstring)
  ref_profile <- Biostrings::DNAStringSet(cleaned_chars)

  # 5. Align ASVs to Profile
  message("  - Aligning ASVs to SILVA Profile using DECIPHER...")
  # Degap ASVs to treat them as unaligned queries
  asv_clean <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(asv_sequences)))
  names(asv_clean) <- names(asv_sequences)

  # Step 1: Align ASVs to each other (temporarily)
  asv_profile <- DECIPHER::AlignSeqs(asv_clean, verbose = FALSE)
  # Step 2: Align the ASV profile to the Reference profile
  merged_alignment <- DECIPHER::AlignProfiles(pattern = asv_profile, subject = ref_profile)

  # Filter to keep only the study sequences (exclude reference profile sequences)
  aligned_asvs <- merged_alignment[names(asv_clean)]

  # 6. Extract Coordinates from Alignment
  message("  - Extracting coordinates...")
  # Helper to find first and last non-gap character
  get_coords <- function(seq_str) {
    chars <- strsplit(seq_str, "")[[1]]
    non_gaps <- which(chars != "-")
    if (length(non_gaps) == 0) return(c(NA, NA))
    return(c(min(non_gaps), max(non_gaps)))
  }

  seq_strs <- as.character(aligned_asvs)
  coords_list <- lapply(seq_strs, get_coords)

  all_ref_starts <- sapply(coords_list, `[`, 1)
  all_ref_ends   <- sapply(coords_list, `[`, 2)

  valid_idx <- !is.na(all_ref_starts) & !is.na(all_ref_ends)
  all_ref_starts <- all_ref_starts[valid_idx]
  all_ref_ends <- all_ref_ends[valid_idx]

  if (length(all_ref_starts) == 0) {
    warning("No valid mapping coordinates found.")
    return(NULL)
  }

  # 7. Calculate Consensus Region for this Study
  consensus <- calculate_consensus_coordinates(all_ref_starts, all_ref_ends, method = config$STUDY_ALIGNMENT_METHOD)
  modal_start <- consensus$start
  modal_end <- consensus$end
  analysis_length <- modal_end - modal_start + 1

  # 8. Build Output Summary
  output_summary <- tibble::tibble(
    study_name = study_info$study_name,
    primer_name = study_info$primer_name,
    ref_start_method = paste0("DECIPHER_SILVA_", config$STUDY_ALIGNMENT_METHOD),
    ref_start = modal_start,
    ref_end = modal_end,
    analysis_amplicon_length = analysis_length,
    original_start_theoretical = NA,
    initial_fwd_trim = study_info$initial_fwd_trim,
    initial_rev_trim = study_info$initial_rev_trim,
    num_asvs_mapped = length(all_ref_starts)
  )

  # 9. Build Detailed Coordinate Table (NEW)
  # Preserves per-ASV data for downstream debugging or detailed analysis
  detailed_coords <- tibble::tibble(
    study_name = study_info$study_name,
    asv_id = names(asv_clean)[valid_idx],
    start = all_ref_starts,
    end = all_ref_ends
  )

  message(paste("  - Done. SILVA Column Region:", modal_start, "-", modal_end))

  return(list(
    summary = output_summary,
    coords = detailed_coords,
    alignment = aligned_asvs
  ))
}

# ------------------------------------------------------------------------------
# Function: relax_consensus_coords
# ------------------------------------------------------------------------------
# Purpose:
#   Safety net that detects and repairs degenerate consensus coordinates where
#   start >= end. This can occur when studies in the overlapping clique have
#   very narrow or contradictory amplicon windows. Without this fix, downstream
#   consumers (visualisation, aggregation, verdict generation) would fail or
#   produce nonsensical results.
#
# Parameters:
#   @param clique_result [list] Output from find_largest_overlapping_clique();
#     must contain $start, $end, and $n_studies.
#   @param all_starts    [integer] Start coordinates for all studies (not just
#     the clique), used as a fallback to define a valid window.
#   @param all_ends      [integer] End coordinates for all studies.
#   @param study_names   [character] Names of all studies (used to update
#     n_studies if relaxation is applied).
#   @param study_name    [character(1)] Label for log messages (default "consensus").
#
# Returns:
#   @return [list] The input clique_result, potentially with $start, $end, and
#     $n_studies corrected. Guarantees start < end with a minimum 50 bp window.
#
# Notes:
#   The 50 bp minimum window matches the min_overlap parameter used elsewhere
#   in the clique-finding algorithm, ensuring consistency.
# ------------------------------------------------------------------------------
relax_consensus_coords <- function(clique_result, all_starts, all_ends, study_names, study_name = "consensus") {
  cons_start <- clique_result$start
  cons_end <- clique_result$end
  
  if (!is.na(cons_start) && !is.na(cons_end) && cons_start >= cons_end) {
    message(sprintf("[WARN] Invalid Consensus Detected (Start: %d >= End: %d) for %s. Applying smart relaxation...", 
                    cons_start, cons_end, study_name))
    
    # Use ALL study bounds as conservative fallback (not just clique)
    study_min <- min(all_starts, na.rm = TRUE)
    study_max <- max(all_ends, na.rm = TRUE)
    
    # Ensure minimum 50bp window matching min_overlap
    relaxed_start <- min(study_min, study_max - 50)
    relaxed_end <- max(study_max, study_min + 50)
    
    clique_result$start <- relaxed_start
    clique_result$end <- relaxed_end
    clique_result$n_studies <- length(study_names)  # Flag as "relaxed consensus"
    
    message(sprintf("  -> Relaxed to study bounds [%d,%d] (length: %d bp)", 
                    relaxed_start, relaxed_end, relaxed_end - relaxed_start))
  }
  return(clique_result)
}
