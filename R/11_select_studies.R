#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   11_select_studies.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    11 - Study Selection
# PURPOSE:  Interactive two-stage CLI for selecting taxonomic level and studies
#
# OVERVIEW:
#   This script provides an interactive wizard that guides the user through
#   study selection for the SeCAT trimming pipeline. First, the user selects
#   a taxonomic level at which to evaluate trim verdicts. Then, verdicts are
#   computed for each study (KEEP, CAUTION_IMPACT, DROP_OUTLIER) based on
#   geometric outlier detection and biological impact thresholds. The user
#   can accept the default (KEEP-only) selection or manually override it.
#   Upon confirmation, the script auto-launches 12_trim_sequences.R.
#
# INPUTS:
#   - output/aggregated_data/verdict_data_all_levels.csv
#       Pre-computed verdict data across all taxonomic levels (from Stage 10)
#
# OUTPUTS:
#   - output/aggregated_data/selected_studies_for_trim.txt
#       Newline-delimited list of study names selected for trimming
#   - output/aggregated_data/selected_analysis_level.txt
#       Single-line file recording the chosen taxonomic level
#   - output/aggregated_data/final_trim_verdicts.csv
#       CSV of study verdicts (Study, Final_Verdict, Reason)
#
# DEPENDENCIES:
#   - dplyr, tidyr, ggplot2, readr, stringr, purrr, tibble, here
#
# CALLED BY:
#   - Manual execution or Nextflow select_studies module
# ==============================================================================

# --- Load Required Packages ---
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

cat("\n================================================================================\n")
cat("                    MESAP FINAL TRIM SELECTION WIZARD                           \n")
cat("================================================================================\n")

# --- 1. Load Raw Verdict Data (All Levels) ---
# The verdict file contains pre-computed trimming assessments for every study
# at every taxonomic level. It is produced by the earlier verdict generation step.
verdict_file <- here("output/aggregated_data/verdict_data_all_levels.csv")

if (!file.exists(verdict_file)) {
  stop("FATAL: Verdict data not found. Run generate_trim_verdicts.R first.")
}

master_verdicts <- read_csv(verdict_file, show_col_types = FALSE)

# --- 2. Stage 1: Select Taxonomic Level ---
# The user chooses which taxonomic rank (e.g., Family, Genus) to use for
# evaluating whether trimming will degrade biological signal. Different ranks
# may yield different verdict outcomes because diversity sensitivity varies.
cat("\n--- STAGE 1: SELECT TAXONOMIC LEVEL ---\n\n")

available_levels <- unique(master_verdicts$Level)
n_levels <- length(available_levels)

cat(sprintf("Found %d taxonomic levels in the analysis:\n\n", n_levels))

# Create a summary table for each level
level_summary <- tibble()
for (i in seq_along(available_levels)) {
  level <- available_levels[i]
  level_data <- master_verdicts %>% filter(Level == !!level)

  unique_studies <- length(unique(level_data$Study))

  level_summary <- bind_rows(level_summary, tibble(
    ID = i,
    Level = level,
    Studies = unique_studies
  ))
}

# Display summary
cat(sprintf("%-3s %-10s %8s\n", "ID", "Level", "Studies"))
cat("-----------------------------\n")

for (i in 1:nrow(level_summary)) {
  row <- level_summary[i,]
  cat(sprintf("%-3d %-10s %8d\n", row$ID, row$Level, row$Studies))
}

cat("\nEnter the Level ID to analyze (e.g., '5' for Family level):\n")
cat("> ")

# Read Level Selection
# Supports both interactive (readline) and non-interactive (stdin) modes
if (interactive()) {
  level_input <- readline()
} else {
  level_input <- readLines("stdin", n=1, warn=FALSE)
}

selected_level_id <- as.integer(trimws(level_input))

if (is.na(selected_level_id) || selected_level_id < 1 || selected_level_id > n_levels) {
  stop(sprintf("Invalid level ID: %s. Must be between 1 and %d.", level_input, n_levels))
}

selected_level <- available_levels[selected_level_id]

cat(sprintf("\nSelected Level: %s\n", selected_level))

# --- 3. Calculate Verdicts for Selected Level ---
# Verdicts are determined by two independent checks per study:
#   (a) Geometric check: Is the study's amplicon region an outlier relative
#       to the consensus? (based on Is_Outlier flag from earlier analysis)
#   (b) Biological impact check: Does the required trimming exceed the
#       changepoint threshold at which diversity loss becomes significant?
# The combined verdict is: DROP_OUTLIER > CAUTION_IMPACT > KEEP.
cat("\nCalculating verdicts for ", selected_level, " level...\n")

level_verdicts <- master_verdicts %>%
  filter(Level == !!selected_level) %>%
  group_by(Study) %>%
  slice(1) %>% # One row per study at this level
  ungroup() %>%
  mutate(
    # 1. Geometric Check
    Status_Geo = if_else(Is_Outlier, "FAIL_OUTLIER", "PASS"),

    # 2. Biological Impact Check
    Safe_Limit = replace_na(Threshold_Observed_Changepoint, Inf),
    Status_Bio = if_else(Threshold_Required >= Safe_Limit, "FAIL_IMPACT", "PASS"),

    # 3. Combined Verdict
    Final_Verdict = case_when(
      Status_Geo == "FAIL_OUTLIER" ~ "DROP_OUTLIER",
      Status_Bio == "FAIL_IMPACT"  ~ "CAUTION_IMPACT",
      TRUE                         ~ "KEEP"
    ),

    # Reason
    Reason = case_when(
      Final_Verdict == "DROP_OUTLIER" ~ "Does not overlap consensus region",
      Final_Verdict == "CAUTION_IMPACT" ~
        paste0("Req. trim (", Threshold_Required, "bp) exceeds safe limit (",
               pmin(Safe_Limit, Threshold_Observed_Changepoint, na.rm=TRUE), "bp)"),
      TRUE ~ "Safe to trim"
    )
  ) %>%
  select(Study, Primer, Final_Verdict, Reason, Threshold_Required, Threshold_Observed_Changepoint) %>%
  arrange(Final_Verdict, Study)

# --- 4. Stage 2: Select Studies (Based on Chosen Level) ---
# Present the user with an ANSI-coloured table of verdicts and allow them
# to confirm the default KEEP-only set or manually specify study IDs.
cat("\n--- STAGE 2: SELECT STUDIES FOR ANALYSIS ---\n\n")

verdicts <- level_verdicts

# Display menu
cat(sprintf("Study Selection at Level: %s\n\n", selected_level))
cat(sprintf("%-4s %-28s %-18s %-45s\n", "ID", "Study Name", "Verdict", "Reason"))
cat("--------------------------------------------------------------------------------\n")

for(i in 1:nrow(verdicts)) {
  row <- verdicts[i,]

  # ANSI Color Codes
  color_start <- switch(row$Final_Verdict,
                       "KEEP" = "\033[32m",           # Green
                       "CAUTION_IMPACT" = "\033[33m", # Yellow
                       "DROP_OUTLIER" = "\033[31m",   # Red
                       "\033[0m")
  color_end <- "\033[0m"

  cat(sprintf("%-4d %-28s %s%-18s%s %-45s\n",
              i,
              substr(row$Study, 1, 28),
              color_start, row$Final_Verdict, color_end,
              substr(row$Reason, 1, 45)))
}
cat("--------------------------------------------------------------------------------\n")

# Get default (KEEP) indices
default_indices <- which(verdicts$Final_Verdict == "KEEP")
default_str <- paste(default_indices, collapse=",")

cat("\n")
cat("Color Legend:\n")
cat("  \033[32m● GREEN (KEEP)\033[0m        - Safe to include\n")
cat("  \033[33m● YELLOW (CAUTION)\033[0m     - May impact diversity; use with caution\n")
cat("  \033[31m● RED (DROP_OUTLIER)\033[0m  - Does not fit consensus\n\n")

n_keep <- sum(verdicts$Final_Verdict == "KEEP")
n_caution <- sum(verdicts$Final_Verdict == "CAUTION_IMPACT")
n_drop <- sum(verdicts$Final_Verdict == "DROP_OUTLIER")

cat(sprintf("Summary: %d KEEP, %d CAUTION, %d DROP_OUTLIER\n\n", n_keep, n_caution, n_drop))

cat("Recommended Action: Proceed with 'KEEP' studies only.\n")
cat(sprintf("Default selection IDs: %s\n", default_str))
cat("Enter IDs to include (comma-separated, e.g. '1,3,5') or press ENTER for default:\n")
cat("> ")

# Read Study Selection
if (interactive()) {
  user_input <- readline()
} else {
  user_input <- readLines("stdin", n=1, warn=FALSE)
}

# --- Parse User Input ---
# Empty input accepts the default KEEP-only set; otherwise parse comma-separated IDs
selected_indices <- integer(0)

if (length(user_input) == 0 || trimws(user_input) == "") {
  # User pressed Enter -> Use defaults
  cat("Using default selection (KEEP only).\n")
  selected_indices <- default_indices
} else {
  # Parse comma-separated string
  clean_input <- gsub(" ", "", user_input)
  parts <- unlist(strsplit(clean_input, ","))

  # Convert to integer
  selected_indices <- as.integer(parts)
  selected_indices <- selected_indices[!is.na(selected_indices)]

  # Validate range
  valid_mask <- selected_indices >= 1 & selected_indices <= nrow(verdicts)
  if (any(!valid_mask)) {
    cat(sprintf("Warning: Ignoring invalid IDs: %s\n", paste(selected_indices[!valid_mask], collapse=",")))
    selected_indices <- selected_indices[valid_mask]
  }
}

# --- 5. Save Selections ---
# Write three output files: the study roster, the taxonomic level choice,
# and the full verdict table for downstream audit/provenance.

if (length(selected_indices) == 0) {
  stop("No valid studies selected. Exiting without saving.")
}

selected_studies <- verdicts$Study[selected_indices]

cat(sprintf("\nSelected %d studies for final trimming:\n", length(selected_studies)))
for (s in selected_studies) {
  cat(sprintf("  - %s\n", s))
}

# Output file paths
studies_file <- here("output/aggregated_data/selected_studies_for_trim.txt")
level_file <- here("output/aggregated_data/selected_analysis_level.txt")
verdicts_file <- here("output/aggregated_data/final_trim_verdicts.csv")

# Write selections
write_lines(selected_studies, studies_file)
write_lines(selected_level, level_file)
write_csv(verdicts %>% select(Study, Final_Verdict, Reason), verdicts_file)

cat(sprintf("\n✓ Studies saved to: %s\n", studies_file))
cat(sprintf("✓ Level saved to:   %s\n", level_file))
cat(sprintf("✓ Verdicts saved to: %s\n", verdicts_file))

# --- 6. Auto-Launch Trimmer ---
# Chain directly into the trimming step (12_trim_sequences.R) so that the
# user experiences a single continuous workflow from selection through trimming.
cat("\n================================================================================\n")
cat("                          LAUNCHING STANDARDIZATION                             \n")
cat("================================================================================\n\n")

trimmer_result <- system("Rscript scripts/12_trim_sequences.R")

if (trimmer_result == 0) {
  cat("\n✅ SELECTION AND STANDARDIZATION COMPLETE!\n")
  cat("================================================================================\n")
} else {
  cat("\n❌ Error during trimming. Check logs above.\n")
  quit(save = "no", status = 1)
}
