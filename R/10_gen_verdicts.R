#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   10_gen_verdicts.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 10 — Final Verdict Export
# PURPOSE:  Export the master verdict table for downstream interactive selection
#
# OVERVIEW:
#   Reads the master_verdict_table.csv produced by 07_aggregate.R and writes
#   it as verdict_data_all_levels.csv. This file contains the KEEP/EXCLUDE
#   decision for every study x taxonomic level combination and serves as the
#   primary input for any downstream interactive study-selection tool or
#   meta-analysis filtering step. The script is deliberately minimal — all
#   statistical logic resides in 07_aggregate.R.
#
# INPUTS:
#   - output/aggregated_data/master_verdict_table.csv — from Stage 7 aggregation
#
# OUTPUTS:
#   - output/aggregated_data/verdict_data_all_levels.csv — per-study, per-level
#     verdicts with consensus status, threshold positions, and quality flags
#
# DEPENDENCIES:
#   - dplyr, tidyr, ggplot2, readr, stringr, purrr, tibble, here
#
# CALLED BY:
#   - Nextflow verdict generation process (Stage 10)
# ==============================================================================

# --- Load libraries ---
suppressPackageStartupMessages({
  suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
  library(here)
})

cat("--- Generating Verdict Data ---\n")

# --- Load the master verdict table produced by 07_aggregate.R ---
master_verdict_file <- here("output/aggregated_data/master_verdict_table.csv")

if (!file.exists(master_verdict_file)) {
  stop("FATAL: master_verdict_table.csv not found. Run aggregation first.")
}

verdicts <- read_csv(master_verdict_file, show_col_types = FALSE)

cat(sprintf("Loaded %d verdict rows covering %d studies and %d levels.\n",
            nrow(verdicts),
            length(unique(verdicts$Study)),
            length(unique(verdicts$Level))))

# --- Write the verdict data for interactive selection / downstream filtering ---
write_csv(verdicts, here("output/aggregated_data/verdict_data_all_levels.csv"))

cat("✓ Verdict data saved for interactive selection.\n")
