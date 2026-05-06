#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   01_primer_mapping.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 1 — Primer Coordinate Mapping (Primer Mode Only)
# PURPOSE:  Map each primer pair to SILVA reference alignment coordinates.
#
# OVERVIEW:
#   To harmonise amplicon datasets generated with different primer pairs
#   (e.g., 515F/806R targeting V4 vs 341F/785R targeting V3-V4), SeCAT must
#   first determine where each primer's amplicon falls within a common
#   reference coordinate space (SILVA 138). This script aligns every unique
#   primer pair from the manifest to the full SILVA reference and records the
#   start/end positions. These coordinates define the theoretical amplicon
#   boundaries and are used downstream to calculate the consensus region —
#   the maximal 16S window shared across all included studies.
#
# INPUTS:
#   - SECAT_MANIFEST_PATH: TSV with columns primer_name, fwd_seq, rev_seq
#   - REFERENCE_DB_PATH: full SILVA aligned FASTA (NR99)
#
# OUTPUTS:
#   - output/intermediate/primer_coords_phase1_output.csv
#     Columns: primer_name, primer_start, primer_end (SILVA positions)
#
# DEPENDENCIES:
#   - R/secat_config.R: global paths and constants
#   - R/secat_utils.R: general utility functions
#   - R/secat_consensus.R: find_consensus_region() alignment logic
#
# CALLED BY:
#   - modules/local/primer_mapping.nf (PRIMER_MAPPING process)
#   - Only used when analysis_mode = "primer"
# ==============================================================================

#===============================================================================
# SECTION 1: INITIALIZATION
#===============================================================================

# --- Load Libraries & Config ---
# Suppress startup messages for cleaner log files in HPC environments
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Biostrings))

# Load global configuration (paths, constants)
source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R")) 

# Load core utility engines
source("R/secat_utils.R")     # General utilities
source("R/secat_consensus.R") # Contains the primer alignment logic

#===============================================================================
# SECTION 2: VALIDATION & DATA LOADING
#===============================================================================

# --- Main ---
message("--- Phase 1: Finding Primer Coordinates using Reference DB ---")

# Check that the reference database exists before attempting to load it
# This prevents obscure errors inside the alignment function
if (!file.exists(REFERENCE_DB_PATH)) {
  stop("FATAL: The main reference database was not found at the path specified in R: ", REFERENCE_DB_PATH)
}

# Load the master manifest containing the primer definitions
message(paste("Loading manifest from:", SECAT_MANIFEST_PATH))
manifest <- read_tsv(SECAT_MANIFEST_PATH, show_col_types = FALSE)

#===============================================================================
# SECTION 3: PRIMER ALIGNMENT
#===============================================================================

# Find primer coordinates using the full reference database
# This calls `find_consensus_region` from R/secat_consensus.R, which:
# 1. Loads the reference DB
# 2. Aligns each unique primer sequence (allowing mismatches)
# 3. Returns the modal start/end positions
primer_coords <- find_consensus_region(manifest, REFERENCE_DB_PATH, VSEARCH_PATH)

# Remove duplicates - keep only unique primer coordinates
# Scientific Note: Multiple studies in the manifest likely use the same
# standard primer sets (e.g., 515F/806R). We only need one coordinate entry
# per unique primer pair.
primer_coords <- primer_coords %>%
  distinct(primer_name, .keep_all = TRUE)

#===============================================================================
# SECTION 4: OUTPUT
#===============================================================================

# Define the output directory and file
output_dir <- "output/intermediate"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_path <- file.path(output_dir, "primer_coords_phase1_output.csv")

# Save the results for downstream use (e.g., Phase 2 simulations)
write_csv(primer_coords, output_path)

message(paste("--- Phase 1 Complete. Primer coordinates saved to:", output_path, "---"))
