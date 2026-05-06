#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   02_study_mapping.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 2 — Individual Study Coordinate Mapping (Worker)
# PURPOSE:  Align one study's ASVs to the SILVA reference and determine its
#           amplicon start/end coordinates in reference alignment space.
#
# OVERVIEW:
#   This is an array-job worker script: each instance processes a single study
#   (identified by task ID) from the manifest. It aligns the study's ASV
#   sequences to the SILVA reference using DECIPHER (profile-to-profile MSA),
#   then determines the study's amplicon boundaries by finding the modal
#   start/end positions across all aligned ASVs. These empirical coordinates
#   replace the theoretical primer-based coordinates (used in primer mode)
#   with actual alignment-derived positions that account for study-specific
#   biases (primer slippage, trimming artefacts, sequencing platform effects).
#
#   The mapping parts from all studies are collected by the COLLECT_MAPS
#   process and merged into study_alignment_coords.csv for consensus
#   region calculation.
#
# INPUTS:
#   - SGE_TASK_ID environment variable (1-indexed study number)
#   - SECAT_MANIFEST_PATH: cleaned manifest TSV
#   - ASV FASTA files referenced in the manifest
#   - REFERENCE_DB_PATH: SILVA aligned FASTA
#
# OUTPUTS:
#   - output/intermediate/study_mapping_parts/mapping_part_<N>.csv
#     (study-level summary: study_name, ref_start, ref_end)
#   - output/intermediate/asv_coordinates/<study>_coords.csv (per-ASV positions)
#   - output/intermediate/aligned_fastas/<study>_aligned.fasta (MSA output)
#
# DEPENDENCIES:
#   - R/secat_config.R: global paths and configuration
#   - R/secat_mapping.R: map_study_to_reference() alignment engine
#   - Biostrings: FASTA I/O
#
# CALLED BY:
#   - modules/local/study_mapping.nf (STUDY_MAPPING process)
# ==============================================================================

# Wrap entire execution in tryCatch for robust error reporting in cluster logs
tryCatch({

  #=============================================================================
  # SECTION 1: INITIALIZATION
  #=============================================================================
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(readr))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(purrr))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(here))
  suppressPackageStartupMessages(library(Biostrings))

  message("=== Loading configuration ===")
  source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R"))    # Global paths
  source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_mapping.R")) # Mapping engine

  #=============================================================================
  # SECTION 2: TASK IDENTIFICATION
  #=============================================================================
  message("=== Getting task ID ===")
  # Retrieve the task ID assigned by the Grid Engine (1-indexed)
  task_id <- as.integer(Sys.getenv("SGE_TASK_ID", unset = NA))
  
  if (is.na(task_id)) {
      stop("FATAL: SGE_TASK_ID must be set. This script is intended to be run via qsub -t 1-N.")
  }

  message(paste("=== Task ID:", task_id, "==="))

  message("=== Loading manifest ===")
  # Load the full manifest to find out which study corresponds to this Task ID
  manifest <- readr::read_tsv(SECAT_MANIFEST_PATH, show_col_types = FALSE)
  
  # Validation: Ensure Task ID is within bounds
  if (task_id > nrow(manifest)) {
      stop(paste("FATAL: Task ID", task_id, "exceeds manifest rows", nrow(manifest)))
  }
  
  # Extract the specific row for this job
  study_info <- manifest[task_id, ]
  message(paste("=== Processing study:", study_info$study_name, "==="))

  #=============================================================================
  # SECTION 3: MAPPING EXECUTION
  #=============================================================================
  message("=== Building config list ===")
  # Package global variables into a clean list for the function API
  config_list <- list(
    ANALYSIS_MODE = ANALYSIS_MODE,
    USE_ALL_ASVS = USE_ALL_ASVS,
    ASV_SAMPLE_SIZE = ASV_SAMPLE_SIZE,
    STUDY_ALIGNMENT_METHOD = STUDY_ALIGNMENT_METHOD,
    REFERENCE_ALIGNMENT_MODE = REFERENCE_ALIGNMENT_MODE,
    REFERENCE_SUBSET_SIZE = REFERENCE_SUBSET_SIZE
  )

  message("=== Calling map_study_to_reference ===")
  # This is the heavy lifting: Subsamples ASVs -> Aligns to SILVA -> Extracts Coords
  mapping_result <- map_study_to_reference(study_info, REFERENCE_DB_PATH, config = config_list)

  #=============================================================================
  # SECTION 4: OUTPUT SAVING
  #=============================================================================
  message("=== Saving output ===")

  # Validate the result structure returned by the engine
  if (is.null(mapping_result)) {
    stop(paste("FATAL: map_study_to_reference returned NULL for", study_info$study_name))
  }

  if (!is.list(mapping_result)) {
    stop(paste("FATAL: map_study_to_reference did not return a list for", study_info$study_name))
  }

  # --- 1. Save Summary (Required) ---
  # This contains the study-level start/end coordinates.
  # We save it as a "part" file; Phase 3 will merge all parts.
  if (!is.null(mapping_result$summary) && is.data.frame(mapping_result$summary)) {
    output_dir <- file.path(OUTDIR, "intermediate", "study_mapping_parts")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    summary_path <- file.path(output_dir, paste0("mapping_part_", task_id, ".csv"))
    readr::write_csv(mapping_result$summary, summary_path)
    message(paste("*** SUCCESS: Summary saved to:", summary_path))
  } else {
    warning("WARNING: Summary is NULL or not a data.frame. Skipping summary save.")
  }

  # --- 2. Save Detailed Coordinates (Optional / Study Mode Only) ---
  # Saves the start/end position for every single ASV processed.
  if (!is.null(mapping_result$coords)) {
    if (is.data.frame(mapping_result$coords)) {
      coords_dir <- file.path(OUTDIR, "intermediate", "asv_coordinates")
      if (!dir.exists(coords_dir)) dir.create(coords_dir, recursive = TRUE)

      coords_path <- file.path(coords_dir, paste0(study_info$study_name, "_coords.csv"))
      readr::write_csv(mapping_result$coords, coords_path)
      message(paste("*** SUCCESS: Coordinates saved to:", coords_path))
    } else {
      warning("WARNING: Coords is not a data.frame. Skipping coords save.")
    }
  } else {
    message("*** INFO: No detailed coordinates returned (Legacy/Primer Mode).")
  }

  # --- 3. Save Full Alignment (Optional / Study Mode Only) ---
  # Saves the MSA of the study ASVs against the reference profile.
  if (!is.null(mapping_result$alignment)) {
    align_dir <- file.path(OUTDIR, "intermediate", "aligned_fastas")
    if (!dir.exists(align_dir)) dir.create(align_dir, recursive = TRUE)

    align_path <- file.path(align_dir, paste0(study_info$study_name, "_aligned.fasta"))
    Biostrings::writeXStringSet(mapping_result$alignment, align_path)
    message(paste("*** SUCCESS: Alignment saved to:", align_path))
  } else {
    message("*** INFO: No alignment returned (Legacy/Primer Mode).")
  }

  message("=== Completed Successfully ===")

}, error = function(e) {
  # Error Handler: Prints full trace to log and ensures non-zero exit code
  message("=== ERROR CAUGHT ===")
  message(paste("Error:", e$message))
  message(paste("Call:", paste(deparse(e$call), collapse = "\n")))
  quit(save = "no", status = 1)
})
