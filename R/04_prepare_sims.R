# ==============================================================================
# SCRIPT:   04_prepare_sims.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 4 - Simulation Planning
# PURPOSE:  Generate the simulation task matrix and SILVA reference subset
#
# OVERVIEW:
#   This script prepares all inputs needed before parallel simulation jobs can
#   execute. It (1) creates a memory-efficient random subset of the full SILVA
#   reference database via reservoir sampling, and (2) builds a task CSV that
#   enumerates every (study-or-primer x replicate-seed) combination the
#   simulation workers (05_sim_worker.R) will process. The task CSV also records
#   the maximum number of trim steps each amplicon requires, ensuring simulations
#   mirror the trimming regime applied to real data.
#
# INPUTS:
#   - SECAT_MANIFEST_PATH          TSV manifest listing studies and primer sets
#   - REFERENCE_DB_PATH            Full SILVA reference FASTA (DNA or RNA)
#   - output/intermediate/primer_coords_phase1_output.csv  (primer mode)
#   - output/intermediate/study_alignment_coords.csv       (study mode)
#
# OUTPUTS:
#   - output/intermediate/simulation_tasks.csv             Master job list
#   - output/intermediate/simulation_reference_subset.fasta Ground-truth community DB
#   - output/intermediate/consensusregioninfo.csv          Consensus region summary (study mode)
#   - output/intermediate/primer_max_trim_steps.csv        Per-primer step counts (primer mode)
#
# DEPENDENCIES:
#   - R packages: dplyr, tidyr, ggplot2, readr, stringr, purrr, tibble, here, Biostrings
#   - Sourced: R/secat_config.R, R/secat_consensus.R (study mode only)
#
# CALLED BY:
#   - Nextflow process: prepare_simulations (runs once before simulation scatter)
# ==============================================================================

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(Biostrings))

# --- Logging Utility ---
# Timestamped messages flushed immediately for real-time monitoring in HPC logs
log_and_flush <- function(message) {
  cat(paste(Sys.time(), "|", message, "\n"))
  flush.console()
}

# ------------------------------------------------------------------------------
# Function: sample_fasta_reservoir
# Purpose:  Draw a uniform random sample of sequences from a large FASTA file
#           without loading the entire file into memory.
#
# Parameters:
#   @param filepath  [character] - Path to FASTA file (can be multi-GB SILVA DB)
#   @param n_samples [integer]   - Number of sequences to retain (reservoir size)
#   @param seed      [integer]   - RNG seed for reproducibility
#
# Returns:
#   @return [list] with elements:
#     $seqs           - character vector of sampled sequences (raw, not DNAStringSet)
#     $headers        - character vector of corresponding FASTA headers
#     $total_scanned  - total number of sequences seen in file
#
# Notes:
#   Uses Vitter's reservoir sampling algorithm: the first k items fill the
#   reservoir directly; each subsequent item i replaces a random existing item
#   with probability k/i, yielding a uniform sample in O(N) time and O(k) space.
#   Returns raw character vectors (not Biostrings objects) to avoid premature
#   RNA/DNA validation errors when SILVA uses 'U' instead of 'T'.
# ------------------------------------------------------------------------------

sample_fasta_reservoir <- function(filepath, n_samples = 100, seed = 123) {
    set.seed(seed)

    log_and_flush(sprintf("  -> Starting reservoir sampling: targeting %d sequences", n_samples))

    con <- file(filepath, "r")
    on.exit(close(con))

    reservoir_seqs <- character(n_samples)
    reservoir_headers <- character(n_samples)
    n_filled <- 0
    total_seen <- 0

    current_header <- NULL
    current_seq_buffer <- character()

    # --- Inner function: decide whether to keep the current sequence ---
    process_sequence <- function(header, seq_fragments) {
        # Concatenate multi-line sequence data into a single string
        full_seq <- paste(seq_fragments, collapse="")

        # Reservoir sampling logic:
        #   Phase 1: Fill reservoir slots sequentially (first n_samples items)
        #   Phase 2: For item i > n_samples, accept with probability n_samples/i
        if (total_seen < n_samples) {
            n_filled <<- n_filled + 1
            reservoir_headers[n_filled] <<- header
            reservoir_seqs[n_filled] <<- full_seq
        } else {
            j <- sample.int(total_seen + 1, 1)
            if (j <= n_samples) {
                reservoir_headers[j] <<- header
                reservoir_seqs[j] <<- full_seq
            }
        }
    }

    # Read file in 50k-line chunks to balance I/O and parsing overhead
    chunk_size <- 50000

    while(TRUE) {
        lines <- readLines(con, n = chunk_size)
        if(length(lines) == 0) break

        # Identify FASTA header lines (">...") within this chunk
        header_indices <- grep("^>", lines)

        if (length(header_indices) == 0) {
            # Entire chunk is continuation of the current sequence
            if (!is.null(current_header)) {
                current_seq_buffer <- c(current_seq_buffer, lines)
            }
            next
        }

        # --- Process chunk containing a mix of headers and sequence lines ---
        prev_idx <- 1
        for (idx in header_indices) {
            # Sequence lines before this header belong to the previous entry
            if (idx > 1) {
                if (!is.null(current_header)) {
                    current_seq_buffer <- c(current_seq_buffer, lines[prev_idx:(idx-1)])
                }
            } else if (idx == 1 && !is.null(current_header)) {
                # Header is first line of chunk; flush the previous sequence
                process_sequence(current_header, current_seq_buffer)
                total_seen <- total_seen + 1
            }

            if (idx > 1 && !is.null(current_header)) {
                process_sequence(current_header, current_seq_buffer)
                total_seen <- total_seen + 1
            }

            # Begin accumulating a new sequence entry
            current_header <- sub("^>", "", lines[idx])
            current_seq_buffer <- character()
            prev_idx <- idx + 1
        }

        # Trailing sequence lines after the last header in this chunk
        if (prev_idx <= length(lines)) {
            current_seq_buffer <- c(current_seq_buffer, lines[prev_idx:length(lines)])
        }

        # Progress indicator (sparse to avoid I/O bottleneck)
        if (total_seen %% 50000 == 0) {
            cat(sprintf("\r    Scanning sequence %d...", total_seen))
        }
    }

    # Flush the final sequence in the file
    if (!is.null(current_header)) {
        process_sequence(current_header, current_seq_buffer)
        total_seen <- total_seen + 1
    }
    cat("\n")

    # Filter out any unfilled reservoir slots (occurs when file has < n_samples seqs)
    valid_indices <- which(reservoir_headers != "")

    log_and_flush(sprintf("  -> Sampling complete: %d sequences selected from %d total",
                          length(valid_indices), total_seen))

    return(list(
        seqs = reservoir_seqs[valid_indices],
        headers = reservoir_headers[valid_indices],
        total_scanned = total_seen
    ))
}

# ------------------------------------------------------------------------------
# Main Logic
# ------------------------------------------------------------------------------
main <- function() {
    log_and_flush("--- SCRIPT STARTED: prepare_simulation_tasks.R V13 ---")

    source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R"))
    log_and_flush(paste("ANALYSIS_MODE detected:", ANALYSIS_MODE))

    # --- Load pipeline configuration with safe defaults ---
    TRIM_INCREMENT <- if (exists("TRIM_INCREMENT")) TRIM_INCREMENT else 10
    MIN_FINAL_SEQUENCE_LENGTH <- if (exists("MIN_FINAL_SEQUENCE_LENGTH")) MIN_FINAL_SEQUENCE_LENGTH else 100
    MAX_ABSOLUTE_TRIM_STEPS <- if (exists("MAX_ABSOLUTE_TRIM_STEPS")) MAX_ABSOLUTE_TRIM_STEPS else 2000
    NUM_SIMULATIONS_PER_PRIMER <- if (exists("NUM_SIMULATIONS_PER_PRIMER")) NUM_SIMULATIONS_PER_PRIMER else 100
    STEP_MODE <- if(exists("TRIM_STEP_MODE")) TRIM_STEP_MODE else "scaled"
    DEFAULT_STEPS <- if(exists("DEFAULT_MAX_TRIM_STEPS")) DEFAULT_MAX_TRIM_STEPS else 40
    BUFFER_STEPS <- if(exists("CONSENSUS_BUFFER_STEPS")) CONSENSUS_BUFFER_STEPS else 10

    # Number of SILVA sequences to include in the simulation ground-truth subset
    SUBSET_SIZE <- if(exists("SIMULATION_MAX_SILVA_SUBSET")) SIMULATION_MAX_SILVA_SUBSET else 10000

    manifest <- readr::read_tsv(SECAT_MANIFEST_PATH, show_col_types = FALSE)
    simulation_tasks <- tibble()

    # =========================================================
    # STEP 1: Generate Reference Subset (RESERVOIR SAMPLING)
    # =========================================================
    # The simulation needs a "ground truth" community drawn from SILVA. Rather
    # than loading the full multi-GB database into every worker, we create a
    # single shared subset here. Each simulation worker then draws its own
    # sub-community from this subset (see 05_sim_worker.R).

    if (TRUE) {
        log_and_flush("\n=== GENERATING SILVA SUBSET FOR SIMULATIONS ===")
        subset_path <- file.path(OUTDIR, "intermediate/simulation_reference_subset.fasta")

        if (!file.exists(subset_path)) {
            log_and_flush(sprintf("  -> Performing Reservoir Sampling (%d seqs) from Reference DB...", SUBSET_SIZE))
            if (!file.exists(REFERENCE_DB_PATH)) stop("FATAL: Reference DB not found at: ", REFERENCE_DB_PATH)

            tryCatch({
                # Draw random subset from full SILVA database
                sample_res <- sample_fasta_reservoir(REFERENCE_DB_PATH, n_samples = SUBSET_SIZE, seed = 42)

                # Convert RNA notation (U) to DNA (T) for VSEARCH compatibility.
                # SILVA stores sequences with RNA alphabet; VSEARCH requires DNA.
                # Also replace non-standard gap characters with standard gaps.
                raw_seqs <- sample_res$seqs
                cleaned_seqs <- gsub(".", "-", raw_seqs, fixed = TRUE)
                cleaned_seqs <- gsub("U", "T", cleaned_seqs, ignore.case = TRUE)

                # Create Biostrings object now that alphabet is safe
                ref_clean <- Biostrings::DNAStringSet(cleaned_seqs)
                names(ref_clean) <- sample_res$headers

                Biostrings::writeXStringSet(ref_clean, subset_path)

                # Report diagnostics
                file_size_mb <- file.info(subset_path)$size / 1024^2
                avg_length <- mean(Biostrings::width(ref_clean))

                log_and_flush(paste("  -> Saved simulation reference subset to:", subset_path))
                log_and_flush(sprintf("  -> Subset contains %d sequences (%.1f MB)",
                                     length(ref_clean), file_size_mb))
                log_and_flush(sprintf("  -> Average sequence length: %.0f bp", avg_length))

            }, error = function(e) {
                stop(paste("Failed to generate reference subset:", e$message))
            })
        } else {
            # --- Validate pre-existing subset ---
            existing_seqs <- Biostrings::readDNAStringSet(subset_path)
            file_size_mb <- file.info(subset_path)$size / 1024^2

            log_and_flush("  -> Simulation reference subset already exists. Skipping generation.")
            log_and_flush(sprintf("  -> Using existing subset: %d sequences (%.1f MB)",
                                 length(existing_seqs), file_size_mb))

            # Warn if existing subset is substantially smaller than configured size
            if (length(existing_seqs) < SUBSET_SIZE * 0.9) {
                warning(sprintf("Existing subset has only %d sequences but config specifies %d. Consider regenerating.",
                               length(existing_seqs), SUBSET_SIZE))
            }
        }
    }

    # =========================================================
    # STEP 2: Task Generation Logic
    # =========================================================
    # Build the simulation task matrix: one row per (task_id x seed) combination.
    # In "primer" mode, task_id = primer_name; in "study" mode, task_id = study_name.
    # Each row tells 05_sim_worker.R the amplicon length and number of trim steps
    # to simulate, ensuring the null model mirrors the real-data trimming regime.

    if (ANALYSIS_MODE == "primer") {
        log_and_flush("\n=== Executing in 'primer' mode. ===")
        coords_path <- file.path(OUTDIR, "intermediate/primer_coords_phase1_output.csv")
        primer_coords <- readr::read_csv(coords_path, show_col_types = FALSE)
        all_primers <- unique(primer_coords$primer_name)

        # --- Helper: compute max trim steps in absolute mode ---
        # Reads actual FASTA files to determine amplicon length, then calculates
        # how many TRIM_INCREMENT-sized steps fit before hitting MIN_FINAL_SEQUENCE_LENGTH.
        # Uses base R subsetting to avoid dplyr non-standard evaluation issues.
        calc_abs <- function(p_name_val) {
             # Explicit subsetting - safe against NSE issues
             studies <- manifest[manifest$primer_name == p_name_val, ]

             if(nrow(studies)==0) return(30)
             max_len <- 0
             for(i in 1:nrow(studies)) {
                 fp <- studies$asv_fasta_path[i]
                 if(file.exists(fp)) {
                     try({ s <- readDNAStringSet(fp); max_len <- max(max_len, max(width(s))) }, silent=TRUE)
                 }
             }
             if(max_len==0) return(30)
             # Max steps = (amplicon_length - minimum_viable_length) / increment
             return(min(floor(max(0, max_len - MIN_FINAL_SEQUENCE_LENGTH)/TRIM_INCREMENT), MAX_ABSOLUTE_TRIM_STEPS))
        }

        # --- Calculate steps per primer ---
        # Scaled mode: fixed number of steps + buffer (increment varies per study)
        # Absolute mode: steps derived from amplicon length and fixed bp increment
        primer_steps_table <- purrr::map_dfr(all_primers, function(curr_primer) {
            steps <- if (STEP_MODE == "scaled") DEFAULT_STEPS + BUFFER_STEPS else calc_abs(curr_primer)
            tibble(primer_name = curr_primer, max_trim_steps = steps)
        })

        readr::write_csv(primer_steps_table, file.path(OUTDIR, "intermediate/primer_max_trim_steps.csv"))

        # --- Build full task matrix using base R ---
        # Avoids dplyr joins that can trigger "object not found" errors in
        # non-interactive Nextflow execution contexts.

        # 1. Compute amplicon length from primer coordinates
        primer_coords$amplicon_length <- primer_coords$primer_end - primer_coords$primer_start + 1
        unique_lengths <- unique(primer_coords[, c("primer_name", "amplicon_length")])

        # 2. Merge step counts with amplicon lengths
        merged_tasks <- merge(primer_steps_table, unique_lengths, by="primer_name")
        colnames(merged_tasks)[colnames(merged_tasks) == "primer_name"] <- "task_id"
        colnames(merged_tasks)[colnames(merged_tasks) == "max_trim_steps"] <- "num_steps"

        # 3. Expand: each task_id x each replicate seed
        seeds <- 1:NUM_SIMULATIONS_PER_PRIMER
        simulation_tasks <- expand.grid(task_id = merged_tasks$task_id, simulation_seed = seeds, stringsAsFactors = FALSE)

        # 4. Re-attach amplicon metadata
        simulation_tasks <- merge(simulation_tasks, merged_tasks, by="task_id")
        simulation_tasks <- simulation_tasks[, c("task_id", "num_steps", "amplicon_length", "simulation_seed")]
        simulation_tasks <- as_tibble(simulation_tasks)

    } else if (ANALYSIS_MODE == "study") {
        log_and_flush("\n=== Executing in 'study' mode. ===")
        coords_path <- file.path(OUTDIR, "intermediate/study_alignment_coords.csv")
        study_coords <- readr::read_csv(coords_path, show_col_types = FALSE)

        # --- Build task matrix: studies x replicates ---
        # Each study gets num_steps determined by its distance from the consensus
        # region. In scaled mode this is a fixed step count; in absolute mode
        # it depends on the study's actual amplicon length.
        simulation_tasks <- study_coords %>%
            dplyr::mutate(
                task_id = study_name,
                num_steps = if (STEP_MODE == "scaled") {
                     DEFAULT_STEPS + BUFFER_STEPS
                } else {
                    # Absolute mode: steps = (amplicon - min_length) / increment
                    steps_calc <- floor((analysis_amplicon_length - MIN_FINAL_SEQUENCE_LENGTH) / TRIM_INCREMENT)
                    pmin(pmax(0, steps_calc), MAX_ABSOLUTE_TRIM_STEPS)
                }
            ) %>%
            dplyr::select(task_id, num_steps, amplicon_length = analysis_amplicon_length) %>%
            tidyr::crossing(simulation_seed = 1:NUM_SIMULATIONS_PER_PRIMER)

        # --- Calculate and save consensus region info ---
        # The consensus region is the maximal 16S region shared across all studies
        # after SILVA alignment. Studies whose amplicons fall outside this region
        # are flagged as outliers. The consensus drives where trimming starts.
        log_and_flush("--- Calculating and saving consensus region info ---")

        # Load the clique-based consensus finder
        source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_consensus.R"))

        # Find the largest set of mutually overlapping studies (clique algorithm)
        # and define the consensus as the intersection of their aligned coordinates
        consres <- find_largest_overlapping_clique(
          starts = study_coords$ref_start,
          ends = study_coords$ref_end,
          study_names = study_coords$study_name,
          min_overlap = 50
        )

        # Summarise consensus region for downstream scripts
        consensus_info <- tibble::tibble(
          ConsensusStart = if(!is.null(consres$start)) consres$start else NA,
          ConsensusEnd = if(!is.null(consres$end)) consres$end else NA,
          ConsensusLength = if(!is.null(consres$start)) (consres$end - consres$start) else 0,
          NumStudiesInConsensus = if(!is.null(consres$n_studies)) consres$n_studies else 0,
          IncludedStudies = if(!is.null(consres$included_studies)) paste(consres$included_studies, collapse = ";") else "",
          OutlierStudies = if(!is.null(consres$excluded_studies)) paste(consres$excluded_studies, collapse = ";") else ""
        )

        # Write consensus info for use by simulation workers and real-data analysis
        consensus_path <- file.path(OUTDIR, "intermediate/consensusregioninfo.csv")
        readr::write_csv(consensus_info, consensus_path)
        log_and_flush(sprintf("✓ Saved consensus region info to %s", consensus_path))

        # Log consensus details for pipeline traceability
        if (!is.null(consres$start)) {
            log_and_flush(sprintf("  Consensus: %d-%d bp (%d studies)",
                                 consres$start, consres$end, consres$n_studies))
            log_and_flush(sprintf("  Outliers: %s",
                                 ifelse(length(consres$excluded_studies) > 0,
                                        paste(consres$excluded_studies, collapse = ", "),
                                        "None")))
        } else {
            log_and_flush("  WARNING: No valid consensus clique could be found.")
        }
    }

    # --- Save Master Task List ---
    # This CSV is the critical handoff to the Nextflow scatter step: each row
    # becomes one parallel invocation of 05_sim_worker.R.
    output_path <- file.path(OUTDIR, "intermediate/simulation_tasks.csv")
    readr::write_csv(simulation_tasks, output_path)
    log_and_flush(sprintf("Saved %d simulation tasks to %s", nrow(simulation_tasks), output_path))
    log_and_flush(sprintf("  -> %d unique tasks × %d replicates each",
                         length(unique(simulation_tasks$task_id)),
                         NUM_SIMULATIONS_PER_PRIMER))
}
# Run Main
tryCatch(main(), error = function(e) {
    log_and_flush(paste("ERROR:", conditionMessage(e)))
    quit(save = "no", status = 1)
})
