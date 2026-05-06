#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   12_trim_sequences.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    12 - Sequence Trimming (Standardisation)
# PURPOSE:  Extract the consensus 16S region from SILVA-aligned ASV sequences
#
# OVERVIEW:
#   This script trims each study's SILVA-aligned FASTA to the global consensus
#   region, then removes alignment gaps to produce degapped sequences that are
#   positionally equivalent across all studies. Trimming is performed in SILVA
#   alignment coordinate space (not raw sequence space), ensuring that column
#   indices correspond to homologous positions regardless of primer pair.
#   Studies that do not overlap the consensus region, or that yield too few
#   sequences above the minimum length threshold, are excluded.
#
# INPUTS:
#   - output/aggregated_data/selected_studies_for_trim.txt
#       Study roster (or auto/roster mode via env vars)
#   - output/intermediate/study_alignment_coords.csv
#       Per-study SILVA alignment start/end coordinates
#   - output/intermediate/consensusregioninfo.csv
#       Global consensus region boundaries (ConsensusStart, ConsensusEnd)
#   - output/intermediate/aligned_fastas/{study}_aligned.fasta
#       SILVA-aligned FASTA files per study
#
# OUTPUTS:
#   - output/standardized_datasets/{study}_standardized.fasta
#       Degapped, trimmed FASTA per study (consensus region only)
#   - output/standardized_datasets/trim_summary.csv
#       Per-study trimming outcome (status, sequence counts, length ranges)
#
# DEPENDENCIES:
#   - dplyr, tidyr, readr, stringr, purrr, tibble, Biostrings
#   - R/secat_config.R (pipeline configuration)
#
# CALLED BY:
#   - Nextflow trim_sequences.nf module (or chained from 11_select_studies.R)
# ==============================================================================

# --- Load Required Packages ---
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(Biostrings)
})

cat("\n=== STARTING STANDARDIZATION (From Aligned Sequences) ===\n")

# --- Load Pipeline Configuration ---
# Configuration is resolved via environment variables (set by Nextflow) with
# fallbacks to secat_config.R defaults. This allows the same script to run
# both inside and outside of the Nextflow workflow.
PROJECTDIR <- Sys.getenv("SECAT_PROJECTDIR", getwd())
config_file <- file.path(PROJECTDIR, "R", "secat_config.R")
if (file.exists(config_file)) {
  source(config_file)
  cat(sprintf("  ✓ Config loaded from: %s\n", config_file))
} else {
  cat("  ⚠️ Config file not found, using defaults\n")
  MIN_REQUIRED_LENGTH <- 50
  OUTDIR <- "."
}

# Path overrides from env vars (set by Nextflow module)
OUTDIR      <- Sys.getenv("SECAT_OUTDIR", OUTDIR)
OUTPUT_DIR  <- file.path(OUTDIR, "standardized_datasets")
ALIGNED_DIR <- file.path(OUTDIR, "intermediate", "aligned_fastas")

# Minimum degapped length (bp) for a trimmed sequence to be retained
if (!exists("MIN_REQUIRED_LENGTH")) MIN_REQUIRED_LENGTH <- 50
# Minimum fraction of sequences that must pass the length filter per study
if (!exists("MIN_YIELD_RATE"))      MIN_YIELD_RATE      <- 0.50

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  OUTDIR:      %s\n", OUTDIR))
cat(sprintf("  OUTPUT_DIR:  %s\n", OUTPUT_DIR))
cat(sprintf("  ALIGNED_DIR: %s\n", ALIGNED_DIR))

# --- 1. Load Inputs ---
# Determine which studies to process. Three selection modes are supported:
#   "auto"   - read verdict file, keep all KEEP studies automatically
#   "roster" - read a user-provided text file listing study names
#   "file"   - read the selection file written by 11_select_studies.R (default)
selection_file <- file.path(OUTDIR, "aggregated_data", "selected_studies_for_trim.txt")
SELECTION_MODE <- Sys.getenv("SECAT_SELECTION_MODE", "file")

if (SELECTION_MODE == "auto") {
  verdict_file <- file.path(OUTDIR, "aggregated_data", "verdict_data_all_levels.csv")
  if (!file.exists(verdict_file)) stop("FATAL: verdict_data_all_levels.csv not found for auto-selection: ", verdict_file)
  verdicts         <- read_csv(verdict_file, show_col_types = FALSE)
  selected_studies <- verdicts %>%
    filter(Final_Verdict == "KEEP") %>%
    pull(Study) %>%
    unique()
  cat(sprintf("Auto-selected %d KEEP studies\n", length(selected_studies)))

} else if (SELECTION_MODE == "roster") {
  # Roster mode: external file path specified via SECAT_SELECTION_FILE env var
  roster_path <- Sys.getenv("SECAT_SELECTION_FILE", "")
  if (!file.exists(roster_path)) stop("FATAL: SECAT_SELECTION_FILE not found: ", roster_path)
  selected_studies <- readLines(roster_path)
  # Strip comment lines and blank lines from the roster
  selected_studies <- selected_studies[!grepl("^#", selected_studies) & nchar(trimws(selected_studies)) > 0]
  cat(sprintf("Roster-selected %d studies from: %s\n", length(selected_studies), roster_path))

} else {
  # Default file mode: read the selection written by Stage 11
  if (!file.exists(selection_file)) stop("FATAL: Selection file not found: ", selection_file)
  selected_studies <- readLines(selection_file)
  cat(sprintf("Loaded %d selected studies from file\n", length(selected_studies)))
}

cat(sprintf("Targeting %d selected studies...\n", length(selected_studies)))

# --- Load Alignment Coordinates and Consensus Region ---
# study_alignment_coords.csv maps each study to its SILVA alignment column range.
# consensusregioninfo.csv defines the global consensus start/end columns that
# all studies will be trimmed to.
coords_file <- file.path(OUTDIR, "intermediate", "study_alignment_coords.csv")
if (!file.exists(coords_file)) stop("FATAL: Coordinates file not found: ", coords_file)
coords <- read_csv(coords_file, show_col_types = FALSE)

consensus_file <- file.path(OUTDIR, "intermediate", "consensusregioninfo.csv")
if (!file.exists(consensus_file)) stop("FATAL: Consensus region file not found: ", consensus_file)
consensus_info <- read_csv(consensus_file, show_col_types = FALSE)

consensus_start <- consensus_info$ConsensusStart[1]
consensus_end   <- consensus_info$ConsensusEnd[1]

if (is.na(consensus_start) || is.na(consensus_end)) {
  stop("FATAL: Consensus coordinates missing or invalid.")
}

cat(sprintf("Global Consensus Target: %d - %d (Length: %d bp)\n\n",
            consensus_start, consensus_end, consensus_end - consensus_start + 1))

# --- 2. Process Each Study ---
# For each study: load aligned FASTA, extract the consensus region columns,
# remove gap characters, filter by minimum length, and write output FASTA.
# Each study can fail at multiple checkpoints (no coordinates, no overlap,
# too short after trimming, low yield rate).
trim_summary <- tibble(
  study_name          = character(),
  status              = character(),
  original_seqs       = integer(),
  trimmed_seqs        = integer(),
  aligned_length      = integer(),
  degapped_length_min = integer(),
  degapped_length_max = integer()
)

for (study in selected_studies) {
  cat(sprintf("Processing: %s\n", study))

  # Look up this study's SILVA alignment coordinate range
  study_coords <- coords %>% filter(study_name == !!study)
  if (nrow(study_coords) == 0) {
    cat("  [SKIP] No coordinates found.\n\n")
    next
  }

  aligned_fasta <- file.path(ALIGNED_DIR, paste0(study, "_aligned.fasta"))

  if (!file.exists(aligned_fasta)) {
    cat(sprintf("  [FAIL] Aligned FASTA not found: %s\n\n", aligned_fasta))
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_NO_ALIGNMENT",
      original_seqs = NA_integer_, trimmed_seqs = 0L,
      aligned_length = NA_integer_, degapped_length_min = NA_integer_,
      degapped_length_max = NA_integer_
    ))
    next
  }

  # Load the SILVA-aligned sequences (all same width within a study)
  aligned_seqs   <- readDNAStringSet(aligned_fasta)
  original_count <- length(aligned_seqs)
  cat(sprintf("  -> Loaded %d aligned sequences\n", original_count))

  alignment_widths <- unique(width(aligned_seqs))
  cat(sprintf("  -> Alignment width(s): %s\n", paste(alignment_widths, collapse = ", ")))

  study_ref_start <- study_coords$ref_start[1]
  study_ref_end   <- study_coords$ref_end[1]

  cat(sprintf("  -> Study maps to SILVA: %d - %d\n", study_ref_start, study_ref_end))
  cat(sprintf("  -> Consensus region:     %d - %d\n", consensus_start, consensus_end))

  # Check whether this study's amplicon region overlaps the consensus at all
  if (study_ref_end < consensus_start || study_ref_start > consensus_end) {
    cat("  [FAIL] Study does not overlap consensus region.\n\n")
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_NO_OVERLAP",
      original_seqs = original_count, trimmed_seqs = 0L,
      aligned_length = alignment_widths[1], degapped_length_min = NA_integer_,
      degapped_length_max = NA_integer_
    ))
    next
  }

  # Compute the effective trim window: the intersection of the study's
  # alignment range and the global consensus region
  effective_start <- max(consensus_start, study_ref_start)
  effective_end   <- min(consensus_end,   study_ref_end)
  aligned_len     <- effective_end - effective_start + 1L

  cat(sprintf("  -> Extracting alignment columns: %d - %d\n", effective_start, effective_end))

  # Extract the consensus columns from the aligned sequences, then degap
  trimmed_aligned  <- subseq(aligned_seqs, start = effective_start, end = effective_end)
  trimmed_degapped <- DNAStringSet(gsub("-", "", as.character(trimmed_aligned)))

  degapped_lengths <- width(trimmed_degapped)
  min_len <- min(degapped_lengths)
  max_len <- max(degapped_lengths)
  cat(sprintf("  -> Degapped length range: %d - %d bp\n", min_len, max_len))

  # Filter out sequences that are too short after degapping
  valid_seqs    <- trimmed_degapped[degapped_lengths >= MIN_REQUIRED_LENGTH]
  dropped_count <- original_count - length(valid_seqs)
  if (dropped_count > 0) {
    cat(sprintf("  -> Dropped %d sequences (< %d bp after degapping)\n",
                dropped_count, MIN_REQUIRED_LENGTH))
  }

  if (length(valid_seqs) == 0L) {
    cat("  [FAIL] All sequences too short after trimming.\n\n")
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_TOO_SHORT",
      original_seqs = original_count, trimmed_seqs = 0L,
      aligned_length = aligned_len,
      degapped_length_min = min_len, degapped_length_max = max_len
    ))
    next
  }

  # Check yield rate: if too many sequences were lost, flag the study
  yield_rate <- length(valid_seqs) / original_count
  if (yield_rate < MIN_YIELD_RATE) {
    cat(sprintf("  [FAIL] Only %.1f%% of sequences passed length filter (minimum: %.0f%%).\n\n",
                100 * yield_rate, 100 * MIN_YIELD_RATE))
    trim_summary <- bind_rows(trim_summary, tibble(
      study_name = study, status = "FAIL_LOW_YIELD",
      original_seqs = original_count, trimmed_seqs = length(valid_seqs),
      aligned_length = aligned_len,
      degapped_length_min = as.integer(min(width(valid_seqs))),
      degapped_length_max = as.integer(max(width(valid_seqs)))
    ))
    next
  }

  # --- Write Standardised Output ---
  output_fasta <- file.path(OUTPUT_DIR, paste0(study, "_standardized.fasta"))
  writeXStringSet(valid_seqs, output_fasta, width = 80)

  final_lengths <- width(valid_seqs)
  cat(sprintf("  [OK] Saved %d sequences (yield: %.1f%%)\n",
              length(valid_seqs), 100 * yield_rate))
  cat(sprintf("  -> Output: %s\n\n", output_fasta))

  trim_summary <- bind_rows(trim_summary, tibble(
    study_name = study, status = "SUCCESS",
    original_seqs = original_count, trimmed_seqs = length(valid_seqs),
    aligned_length = aligned_len,
    degapped_length_min = as.integer(min(final_lengths)),
    degapped_length_max = as.integer(max(final_lengths))
  ))
}

# --- 3. Save Summary ---
summary_file <- file.path(OUTPUT_DIR, "trim_summary.csv")
write_csv(trim_summary, summary_file)
cat(sprintf("✓ Trimming summary saved: %s\n", summary_file))

# --- 4. Print Summary ---
cat("\n================================================================================\n")
cat("                        TRIMMING SUMMARY\n")
cat("================================================================================\n\n")
print(trim_summary, n = Inf)

n_success   <- sum(trim_summary$status == "SUCCESS",         na.rm = TRUE)
n_low_yield <- sum(trim_summary$status == "FAIL_LOW_YIELD",  na.rm = TRUE)
n_failed    <- sum(startsWith(trim_summary$status, "FAIL"),  na.rm = TRUE) - n_low_yield

cat(sprintf("\n✓ Successfully trimmed : %d / %d studies\n",
            n_success, length(selected_studies)))
if (n_low_yield > 0) {
  low_yield_studies <- trim_summary$study_name[trim_summary$status == "FAIL_LOW_YIELD"]
  cat(sprintf("⚠ Low yield (excluded): %d — %s\n",
              n_low_yield, paste(low_yield_studies, collapse = ", ")))
}
if (n_failed > 0) {
  cat(sprintf("✗ Other failures       : %d\n", n_failed))
}

cat("\n=== STANDARDIZATION COMPLETE ===\n")
# NOTE: Merger is called as a separate Nextflow process (MERGE_DATASETS),
# not via system() from this script.
