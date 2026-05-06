#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   08_render_wrapper.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 8 (auxiliary) — Report Rendering Wrapper
# PURPOSE:  Thin wrapper that invokes rmarkdown::render for a single study
#
# OVERVIEW:
#   This script acts as a robust bridge between the Nextflow/SGE execution
#   environment and RMarkdown. It accepts a study name as a command-line
#   argument and passes it to rmarkdown::render with correctly resolved
#   absolute paths. This indirection avoids shell quoting and parameter-
#   passing issues that arise when calling rmarkdown directly from a
#   scheduler job script.
#
# INPUTS:
#   - Command-line argument: study_name (character string)
#   - scripts/generate_report.Rmd — the RMarkdown template for per-study reports
#
# OUTPUTS:
#   - output/reports/<study_name>_report.html — rendered per-study HTML report
#
# DEPENDENCIES:
#   - rmarkdown, here
#
# CALLED BY:
#   - Nextflow report generation process (Stage 8)
# ==============================================================================

# --- Load path resolution library ---
suppressPackageStartupMessages(library(here))

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("FATAL: No study_name was provided to the render wrapper script.")
}
study_name_param <- args[1]

# --- Log receipt of parameter for traceability ---
message(paste("Wrapper script received study_name:", study_name_param))

# --- Invoke rmarkdown::render with absolute paths ---
# Using here() ensures paths resolve correctly regardless of the working
# directory set by the scheduler or Nextflow executor.
rmarkdown::render(
  input = here("scripts", "generate_report.Rmd"),
  params = list(study_name = study_name_param),
  output_file = here("output", "reports", paste0(study_name_param, "_report.html"))
)

message("Report rendering complete.")

