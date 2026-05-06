#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   validation_taxon_v2.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 11 -- Multi-Tier Ecological Validation
# PURPOSE:  Validate the harmonised (merged) dataset by comparing pre- and
#           post-trimming community structure across multiple ecological metrics.
#
# OVERVIEW:
#   After trimming and merging, this script assesses whether the harmonisation
#   preserved the ecological signal of the original datasets. It runs a
#   multi-tier validation at each requested taxonomic level (ASV, Genus, Family):
#
#   Tier 0:  Sequence-level QC -- length distributions, GC content stability,
#            and consensus coordinate verification against trim_summary.csv.
#   Tier 1A: Alpha diversity preservation -- ICC (two-way agreement), Spearman
#            rank correlation, and Bland-Altman agreement plots for Observed
#            richness, Chao1, ACE, Shannon, Simpson, and Pielou evenness.
#   Tier 1B: Beta diversity preservation -- Mantel test (Pearson on Bray-Curtis
#            distance matrices), Procrustes/PROTEST rotation analysis, and
#            delta-R2 PERMANOVA comparing study-level variance before/after.
#   Tier 1C: Within-group dispersion -- betadisper (Levene analogue for
#            multivariate data) tests whether trimming shifts within-study
#            variance in ordination space.
#   Tier 1D: Rank abundance curve preservation -- per-sample Spearman rho on
#            sorted relative abundances; tests dominance gradient stability.
#   Tier 2A: Cross-study feature sharing -- MetaASV sharing rate at ASV level
#            (via asv_mapping_final.tsv), pairwise Jaccard at Genus/Family.
#   Tier 3A: Taxonomic composition -- stacked bar plots of top-10 taxa with
#            consistent colour palette across before/after panels.
#   Tier 3B: Co-occurrence network stability -- topology metrics (vertices,
#            edges, density, transitivity), hub node retention (top 10% by
#            degree), and optional meconetcomp robustness analysis.
#   Tier 4:  Abundance concordance -- per-taxon Spearman rho across samples
#            (replaces DESeq2). Uses MetaASV bridge at ASV level.
#
#   Results are saved as CSV tables and PDF plots in the validation output
#   directory. These provide evidence for thesis examiners and reviewers that
#   SeCAT's trimming did not introduce systematic bias.
#
# INPUTS:
#   - Post-trimming: feature_table.tsv, taxonomy.tsv, metadata.tsv,
#     sequences.fasta (in post_consensus/)
#   - Pre-trimming: matching tables in pre_consensus/
#   - ASV mapping (asv_mapping_final.tsv -- old-to-new ASV ID correspondence)
#   - Consensus region info and trim summary CSVs from pipeline output
#
# OUTPUTS:
#   - outputs/ directory: validation_summary.csv, per-tier CSVs, per-level
#     subdirectories with detailed results
#   - outputs/figures/ directory: PDF diagnostic plots per tier and level
#   - outputs/checkpoints/ directory: RDS checkpoints for resumable runs
#
# DEPENDENCIES:
#   - vegan:      diversity indices, PERMANOVA (adonis2), Mantel, betadisper,
#                 PCoA (cmdscale), Procrustes/PROTEST
#   - irr:        intraclass correlation coefficient (ICC)
#   - ggplot2:    all diagnostic visualisations
#   - patchwork:  multi-panel plot assembly
#   - microeco:   microtable objects, trans_beta, trans_abund, trans_network
#   - igraph:     co-occurrence network topology metrics
#   - ape:        phylogenetic utilities (cmdscale wrapper)
#   - Biostrings: (optional) sequence length and GC content analysis
#   - meconetcomp:(optional) network comparison, robustness, Venn diagrams
#
# CALLED BY:
#   - modules/local/validate.nf (VALIDATE process in Nextflow pipeline)
#
# USAGE:
#   Rscript validation_taxon_v2.R [START_TIER] [INSTALL_PACKAGES]
#   Rscript validation_taxon_v2.R /path/to/base START_TIER INSTALL_PACKAGES
#
# Author:  SeCAT Development Team
# Version: 4.1 (2026-04)
# ==============================================================================

# ==============================================================================
# ARGUMENT PARSING
# ------------------------------------------------------------------------------
# Supports three calling conventions:
#   0 args: run from current directory, start at Tier 0, auto-install packages
#   1 arg:  specify START_TIER (integer) to resume from a specific validation tier
#   2 args: START_TIER + INSTALL_PACKAGES (logical, FALSE to skip installs)
#   3 args: BASE_DIR + START_TIER + INSTALL_PACKAGES (full control)
# The checkpoint system (see save_checkpoint/load_checkpoint below) allows
# resuming from any tier without re-loading data.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  BASE_DIR         <- getwd()
  START_TIER       <- 0
  INSTALL_PACKAGES <- TRUE
} else if (length(args) == 1) {
  BASE_DIR         <- getwd()
  START_TIER       <- as.integer(args[1])
  INSTALL_PACKAGES <- TRUE
} else if (length(args) == 2) {
  BASE_DIR         <- getwd()
  START_TIER       <- as.integer(args[1])
  INSTALL_PACKAGES <- as.logical(args[2])
} else {
  BASE_DIR         <- args[1]
  START_TIER       <- as.integer(args[2])
  INSTALL_PACKAGES <- as.logical(args[3])
}

cat("\n================================================================================\n")
cat("=== SeCAT VALIDATION v4.0 ===\n")
cat("================================================================================\n")
cat("Date:          ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("Base directory:", BASE_DIR, "\n")
cat("Start tier:    ", START_TIER, "\n")
cat("Install pkgs:  ", INSTALL_PACKAGES, "\n")
cat("================================================================================\n\n")

# ==============================================================================
# PACKAGE MANAGEMENT
# ------------------------------------------------------------------------------
# Checks for all required and optional packages, installing missing ones from
# CRAN or Bioconductor as needed. Uses a user-local library path to avoid
# permission issues in shared/container environments.
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(BioC_mirror = "https://bioconductor.org")

lib_path <- Sys.getenv("R_LIBS_USER", "~/R/library/4.3")
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))

# check_and_install()
# Purpose:    Test whether a package is available; install if missing and allowed.
# Parameters: pkg      -- package name (string)
#             method   -- "CRAN" or "Bioconductor" (install source)
#             required -- if TRUE, missing package causes pipeline failure
# Returns:    TRUE if package is now available, FALSE otherwise
check_and_install <- function(pkg, method = "CRAN", required = TRUE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %-22s [OK]\n", pkg)); return(TRUE)
  }
  if (!INSTALL_PACKAGES) {
    if (required) cat(sprintf("  ✗ %-22s [MISSING]\n", pkg))
    else          cat(sprintf("  ⚠ %-22s [OPTIONAL MISSING]\n", pkg))
    return(!required)
  }
  cat(sprintf("  ⚙ %-22s [INSTALLING...]\n", pkg))
  tryCatch({
    if (method == "CRAN") {
      install.packages(pkg, lib = lib_path, dependencies = TRUE, quiet = FALSE)
    } else {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", lib = lib_path, quiet = FALSE)
      BiocManager::install(pkg, lib = lib_path, update = FALSE, ask = FALSE)
    }
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  ✓ %-22s [INSTALLED]\n", pkg)); return(TRUE)
    } else {
      cat(sprintf("  ✗ %-22s [INSTALL FAILED]\n", pkg)); return(FALSE)
    }
  }, error = function(e) {
    cat(sprintf("  ✗ %-22s [ERROR: %s]\n", pkg, conditionMessage(e))); return(FALSE)
  })
}

# Install/verify all dependencies. Required packages cause fatal error if
# missing; optional packages (Biostrings, meconetcomp) degrade gracefully.
cat("Checking packages:\n")
ok <- TRUE
ok <- check_and_install("tidyverse", required = FALSE)   && ok
ok <- check_and_install("vegan")       && ok
ok <- check_and_install("irr")         && ok
ok <- check_and_install("ggplot2")     && ok
ok <- check_and_install("patchwork")   && ok
ok <- check_and_install("RColorBrewer") && ok
ok <- check_and_install("data.table")  && ok
ok <- check_and_install("magrittr")    && ok
ok <- check_and_install("ape")         && ok
ok <- check_and_install("igraph")      && ok
ok <- check_and_install("microeco")    && ok
check_and_install("Biostrings", method = "Bioconductor", required = FALSE)
check_and_install("meconetcomp", required = FALSE)

if (!ok) stop("Missing required packages. Exiting.")
cat("\n✓ All required packages available\n\n")

suppressPackageStartupMessages({
  suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
  library(vegan)
  library(irr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(data.table)
  library(magrittr)
  library(ape)
  library(igraph)
  library(microeco)
  library(parallel)
})

# Feature flags for optional packages -- enables graceful degradation
# Biostrings: needed for Tier 0A/0B (sequence length and GC content analysis)
# meconetcomp: needed for extended network comparison in Tier 3B
HAS_BIOSTRINGS    <- requireNamespace("Biostrings",   quietly = TRUE)
HAS_MECONETCOMP   <- requireNamespace("meconetcomp",  quietly = TRUE)
if (HAS_BIOSTRINGS)  suppressPackageStartupMessages(library(Biostrings))
if (HAS_MECONETCOMP) suppressPackageStartupMessages(library(meconetcomp))

# ==============================================================================
# PATHS
# ------------------------------------------------------------------------------
# Directory layout expected by this script:
#   BASE_DIR/
#     pre_consensus/          -- original (untrimmed) study data
#       feature_table.tsv     -- ASV count matrix (features x samples)
#       taxonomy.tsv          -- QIIME2-style taxonomy (ASV_ID, Taxon, Confidence)
#       metadata.tsv          -- sample metadata with SampleID and StudyID columns
#       sequences.fasta       -- (optional) representative sequences
#     post_consensus/         -- SeCAT-trimmed and merged data (same file formats)
#     asv_mapping_final.tsv   -- maps original ASV hashes to MetaASV IDs
#     output/intermediate/    -- consensusregioninfo.csv (from alignment stage)
#     output/standardized_datasets/ -- trim_summary.csv (from trimming stage)
#     outputs/                -- validation results written here
#       figures/              -- PDF plots
#       checkpoints/          -- RDS checkpoints for resumable runs
# ==============================================================================

PRE_CONSENSUS_DIR  <- file.path(BASE_DIR, "pre_consensus")
PRE_OTU  <- file.path(PRE_CONSENSUS_DIR, "feature_table.tsv")
PRE_TAX  <- file.path(PRE_CONSENSUS_DIR, "taxonomy.tsv")
PRE_META <- file.path(PRE_CONSENSUS_DIR, "metadata.tsv")
PRE_FA   <- file.path(PRE_CONSENSUS_DIR, "sequences.fasta")   # optional

POST_CONSENSUS_DIR <- file.path(BASE_DIR, "post_consensus")
POST_OTU    <- file.path(POST_CONSENSUS_DIR, "feature_table.tsv")
POST_TAX    <- file.path(POST_CONSENSUS_DIR, "taxonomy.tsv")
POST_META   <- file.path(POST_CONSENSUS_DIR, "metadata.tsv")
POST_FA     <- file.path(POST_CONSENSUS_DIR, "sequences.fasta")     # optional

# SeCAT pipeline outputs for coordinate verification
CONSENSUS_INFO <- file.path(BASE_DIR, "output/intermediate/consensusregioninfo.csv")
TRIM_SUMMARY   <- file.path(BASE_DIR, "output/standardized_datasets/trim_summary.csv")
ASV_MAPPING    <- file.path(BASE_DIR, "asv_mapping_final.tsv")

# Taxonomic resolution levels for the validation loop. Override via environment
# variable SECAT_VALIDATION_LEVELS (comma-separated, e.g., "ASV,Genus").
# Default: ASV, Genus, Family -- covers fine-grained through coarse resolution.
TAXONOMIC_LEVELS <- {
  env_levels <- Sys.getenv("SECAT_VALIDATION_LEVELS", "")
  if (nchar(env_levels) > 0) {
    strsplit(env_levels, ",")[[1]]
  } else {
    c("ASV", "Genus", "Family")
  }
}

OUTPUT_DIR     <- file.path(BASE_DIR, "outputs")
PLOTS_DIR      <- file.path(OUTPUT_DIR, "figures")
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "checkpoints")

dir.create(OUTPUT_DIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,      recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

ncores <- max(1, detectCores() - 2)
cat("CPU cores available:", ncores, "\n\n")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# Data loading helpers
# ------------------------------------------------------------------------------

# load_otu_table()
# Purpose:    Read a QIIME2-style feature table (TSV), handling both
#             feature-by-sample and sample-by-feature orientations.
# Parameters: filepath -- path to the TSV file
# Returns:    data.frame with features as rows, samples as columns, numeric
#             counts. Zero-sum features are removed. Missing/duplicate IDs
#             are cleaned automatically.
# Notes:      Auto-detects orientation by checking whether the first column
#             header matches known feature ID names (e.g., "#OTU ID",
#             "Feature ID"). If not, assumes sample-by-feature and transposes.
load_otu_table <- function(filepath) {
  cat("  Loading OTU table:", filepath, "\n")
  first_line <- readLines(filepath, n = 1)
  otu <- if (grepl("^#", first_line)) {
    read_tsv(filepath, skip = 1, show_col_types = FALSE)
  } else {
    read_tsv(filepath, show_col_types = FALSE)
  }

  id_col    <- colnames(otu)[1]
  id_values <- otu[[1]]

  n_miss <- sum(is.na(id_values) | id_values == "")
  if (n_miss > 0) {
    cat("    Removing", n_miss, "missing IDs\n")
    otu <- otu[!is.na(id_values) & id_values != "", ]
    id_values <- otu[[1]]
  }
  n_dup <- sum(duplicated(id_values))
  if (n_dup > 0) {
    cat("    Making", n_dup, "duplicate IDs unique\n")
    otu[[1]] <- make.unique(as.character(id_values))
  }

  feature_first_cols <- c("OTU ID", "#OTU ID", "Feature ID", "#Feature ID",
                           "ASV_ID", "featureid", "id", "#id")
  mat <- as.matrix(otu %>% select(-1))
  mode(mat) <- "numeric"
  rownames(mat) <- as.character(otu[[1]])

  if (!id_col %in% feature_first_cols) {
    cat("    Detected sample-by-feature format — transposing\n")
    mat <- t(mat)
  }

  mat[is.na(mat)] <- 0
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  cat(sprintf("    %d features × %d samples | total reads: %s\n",
              nrow(mat), ncol(mat), format(sum(mat), big.mark = ",")))
  as.data.frame(mat)
}

# load_taxonomy()
# Purpose:    Read and parse a QIIME2-style taxonomy file, splitting the
#             semicolon-delimited lineage string into Kingdom through Species.
# Parameters: filepath -- path to taxonomy TSV
# Returns:    data.frame with ASV IDs as rownames and columns for each
#             taxonomic rank (Kingdom, Phylum, Class, Order, Family, Genus,
#             Species), plus the original Taxon string and Confidence score.
# Notes:      Strips QIIME2 prefix notation (e.g., "g__Vibrio" -> "Vibrio").
#             Handles multiple input formats (ASV_ID/Taxon, ASV_ID/Meta_ID).
load_taxonomy <- function(filepath) {
  cat("  Loading taxonomy:", filepath, "\n")
  tax <- read_tsv(filepath, show_col_types = FALSE)

  # Detect and standardise format
  if (all(c("ASV_ID", "Taxon") %in% colnames(tax))) {
    tax_std <- tax %>% select(ASV_ID, Taxon, any_of("Confidence"))
  } else if (all(c("ASV_ID", "Meta_ID") %in% colnames(tax))) {
    tax_std <- data.frame(ASV_ID = tax$Meta_ID,
                          Taxon  = tax$ASV_ID,
                          Confidence = tax$Confidence,
                          stringsAsFactors = FALSE)
  } else {
    colnames(tax)[1] <- "ASV_ID"
    colnames(tax)[2] <- "Taxon"
    tax_std <- tax %>% select(ASV_ID, Taxon, any_of("Confidence"))
  }

  tax_split <- tax_std %>%
    separate(Taxon,
             into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
             sep = ";", fill = "right", remove = FALSE) %>%
    mutate(across(Kingdom:Species,
                  ~str_trim(str_remove(.x, "^[kpcofgsd]__")))) %>%
    column_to_rownames("ASV_ID")

  cat(sprintf("    %d ASVs | %d unique genera\n",
              nrow(tax_split),
              length(unique(tax_split$Genus[!is.na(tax_split$Genus) & tax_split$Genus != ""]))))
  tax_split
}

# load_metadata()
# Purpose:    Read sample metadata, standardising the sample ID column name.
# Parameters: filepath -- path to metadata TSV
# Returns:    data.frame with sample IDs as rownames. The StudyID column
#             (if present) is used throughout for study-level analyses.
# Notes:      Recognises common QIIME2 sample ID synonyms (#SampleID,
#             sample-id, sampleid, etc.) and renames to "SampleID".
load_metadata <- function(filepath) {
  cat("  Loading metadata:", filepath, "\n")
  meta <- read_tsv(filepath, show_col_types = FALSE)
  id_synonyms <- c("sample-id", "#SampleID", "sampleid", "SampleID",
                   "sample_id", "#sample-id")
  for (syn in id_synonyms) {
    if (syn %in% colnames(meta)) { colnames(meta)[colnames(meta) == syn] <- "SampleID"; break }
  }
  if (!"SampleID" %in% colnames(meta)) colnames(meta)[1] <- "SampleID"
  meta <- meta %>% column_to_rownames("SampleID")
  cat(sprintf("    %d samples | %d variables\n", nrow(meta), ncol(meta)))
  meta
}

# ------------------------------------------------------------------------------
# Taxonomic aggregation
# ------------------------------------------------------------------------------

# aggregate_to_level()
# Purpose:    Collapse ASV-level counts to a higher taxonomic rank (e.g., Genus
#             or Family) by summing counts across all ASVs sharing the same
#             classification at that rank.
# Parameters: dataset_asv -- microeco microtable object at ASV resolution
#             taxa_level  -- target rank ("ASV", "Genus", or "Family")
# Returns:    New microtable object at the requested resolution. If taxa_level
#             is "ASV", returns a deep clone of the input (no aggregation).
# Notes:      Unclassified/unassigned ASVs are excluded from the aggregated
#             table. Uses base R rowsum() for efficient grouped summation.
aggregate_to_level <- function(dataset_asv, taxa_level) {
  if (taxa_level == "ASV") return(dataset_asv$clone(deep = TRUE))
  cat(sprintf("    Aggregating to %s level...\n", taxa_level))

  otu <- dataset_asv$otu_table
  tax <- dataset_asv$tax_table

  if (!taxa_level %in% colnames(tax))
    stop(sprintf("Column '%s' not found in taxonomy table", taxa_level))

  taxa_names <- tax[[taxa_level]]
  valid      <- !is.na(taxa_names) & taxa_names != "" & taxa_names != "unassigned"
  agg_counts <- rowsum(otu[valid, , drop = FALSE], taxa_names[valid], na.rm = TRUE)

  tax_agg <- data.frame(Taxa = rownames(agg_counts),
                        row.names = rownames(agg_counts),
                        stringsAsFactors = FALSE)
  colnames(tax_agg)[1] <- taxa_level

  ds <- microtable$new(otu_table    = as.data.frame(agg_counts),
                       tax_table    = tax_agg,
                       sample_table = dataset_asv$sample_table)
  ds$tidy_dataset()
  cat(sprintf("    → %d %s × %d samples\n",
              nrow(ds$otu_table), taxa_level, ncol(ds$otu_table)))
  ds
}

# make_taxa_unique()
# Purpose:    Ensure no duplicate taxa names at the current resolution level.
#             Duplicates can arise from ambiguous classifications.
# Parameters: dataset    -- microtable object
#             taxa_level -- column to check for duplicates
# Returns:    Modified microtable (in-place update of tax_table)
make_taxa_unique <- function(dataset, taxa_level) {
  tax <- dataset$tax_table
  if (taxa_level %in% colnames(tax)) {
    nm <- tax[[taxa_level]]
    if (any(duplicated(nm))) {
      tax[[taxa_level]] <- make.unique(as.character(nm), sep = "_")
      dataset$tax_table <- tax
    }
  }
  dataset
}

# ------------------------------------------------------------------------------
# Bland-Altman agreement plot
# ------------------------------------------------------------------------------

# bland_altman_plot()
# Purpose:    Generate a Bland-Altman plot for assessing measurement agreement
#             between pre- and post-trimming diversity metrics.
# Parameters: x     -- numeric vector of "before" values
#             y     -- numeric vector of "after" values
#             title -- plot title
#             xlab, ylab -- axis labels
# Returns:    ggplot object showing difference vs mean, with bias line (blue)
#             and 95% limits of agreement (red dashed, mean +/- 1.96*SD).
# Context:    In clinical method comparison (Bland & Altman, 1986), the LoA
#             define the range within which 95% of differences are expected.
#             For ecological validation, points within the LoA indicate that
#             trimming did not introduce systematic bias for that metric.
#             Ideal result: bias near zero, narrow LoA.
bland_altman_plot <- function(x, y, title = "Bland-Altman",
                              xlab = "Mean", ylab = "Difference (After - Before)") {
  means <- (x + y) / 2
  diffs <- y - x
  md    <- mean(diffs, na.rm = TRUE)
  sd_d  <- sd(diffs,  na.rm = TRUE)
  loa_u <- md + 1.96 * sd_d
  loa_l <- md - 1.96 * sd_d

  ggplot(data.frame(m = means, d = diffs), aes(x = m, y = d)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_hline(yintercept = md,    color = "blue",  linetype = "solid",  linewidth = 0.9) +
    geom_hline(yintercept = loa_u, color = "red",   linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = loa_l, color = "red",   linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = 0,     color = "grey50", linetype = "dotted") +
    annotate("text", x = Inf, y = md,    label = sprintf("Bias: %.3f",   md),    hjust=1.1, vjust=-0.5, color="blue", size=3.5) +
    annotate("text", x = Inf, y = loa_u, label = sprintf("+1.96SD: %.3f", loa_u), hjust=1.1, vjust=-0.5, color="red",  size=3.5) +
    annotate("text", x = Inf, y = loa_l, label = sprintf("-1.96SD: %.3f", loa_l), hjust=1.1, vjust=1.5,  color="red",  size=3.5) +
    labs(title = title, x = xlab, y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# ------------------------------------------------------------------------------
# Checkpointing (resumable execution)
# ------------------------------------------------------------------------------

# save_checkpoint() / load_checkpoint()
# Purpose:    Persist intermediate state to RDS files so the validation can
#             resume from any tier without re-loading and re-processing data.
#             Particularly useful when network analysis (Tier 3B) is slow.
# Parameters: id        -- checkpoint identifier string (e.g., "data")
#             data_list -- named list of objects to persist
# Returns:    load_checkpoint returns the stored list, or NULL if not found
save_checkpoint <- function(id, data_list) {
  saveRDS(data_list, file.path(CHECKPOINT_DIR, sprintf("ckpt_%s.rds", id)))
  cat(sprintf("  ✓ Checkpoint saved [%s]\n", id))
}

load_checkpoint <- function(id) {
  fp <- file.path(CHECKPOINT_DIR, sprintf("ckpt_%s.rds", id))
  if (file.exists(fp)) { cat(sprintf("  ↩ Loading checkpoint [%s]\n", id)); readRDS(fp) }
  else NULL
}

# ------------------------------------------------------------------------------
# Summary results accumulator
# ------------------------------------------------------------------------------
# Global tibble that collects pass/fail results from every tier and level.
# Each row records: taxonomic Level, Tier ID, Analysis name, Metric name,
# numeric Value, pass Criteria description, and Status (PASS/MODERATE).
# Written to validation_summary.csv at the end for thesis reporting.

summary_results <- tibble(
  Level    = character(),
  Tier     = character(),
  Analysis = character(),
  Metric   = character(),
  Value    = numeric(),
  Criteria = character(),
  Status   = character()
)

# add_result()
# Purpose:    Append a single validation result to the global summary_results.
# Parameters: level     -- taxonomic level (e.g., "ASV", "Genus", "ALL")
#             tier      -- tier identifier (e.g., "1A", "1B")
#             analysis  -- analysis category name
#             metric    -- specific metric name
#             value     -- numeric result
#             criteria  -- human-readable pass criterion description
#             pass_flag -- logical, TRUE = PASS, FALSE = MODERATE
# Side effect: Modifies summary_results in the global environment (<<-)
add_result <- function(level, tier, analysis, metric, value, criteria, pass_flag) {
  status <- if (pass_flag) "PASS" else "MODERATE"
  summary_results <<- bind_rows(summary_results, tibble(
    Level = level, Tier = tier, Analysis = analysis,
    Metric = metric, Value = round(value, 4),
    Criteria = criteria, Status = status
  ))
}

# ==============================================================================
# TIER 0: SEQUENCE-LEVEL QC + DATA LOADING
# ==============================================================================
# This block runs once (not per taxonomic level). It performs three sequence-
# level checks before loading the count data:
#
#   Tier 0A: Length distributions -- compares the coefficient of variation (CV)
#            of sequence lengths before and after trimming. CV should decrease
#            post-trim because SeCAT aligns sequences to a consensus region,
#            producing more uniform lengths. An increase in CV would suggest
#            the trimming introduced length artefacts.
#
#   Tier 0B: GC content preservation -- GC proportion should remain stable
#            (< 2 percentage-point shift). A systematic GC shift would indicate
#            that trimming preferentially removed AT-rich or GC-rich regions,
#            potentially biasing downstream taxonomic classification.
#
#   Tier 0C: Consensus coordinate verification -- checks that the trim
#            coordinates from the pipeline produced sequences within +/-10 bp
#            of the expected consensus region length, and that >= 90% of
#            studies were trimmed successfully.
# ==============================================================================

if (START_TIER == 0) {

  cat("================================================================================\n")
  cat("STEP 0: SEQUENCE-LEVEL QC\n")
  cat("================================================================================\n\n")

  seq_qc_results <- list()

  # --- 0A: Length distributions ---
  # Requires Biostrings package and FASTA files in both pre/post directories
  if (HAS_BIOSTRINGS && file.exists(PRE_FA) && file.exists(POST_FA)) {
    cat("Tier 0A: Length distributions\n")
    seqs_before <- Biostrings::readDNAStringSet(PRE_FA)
    seqs_after  <- Biostrings::readDNAStringSet(POST_FA)

    len_before <- Biostrings::width(seqs_before)
    len_after  <- Biostrings::width(seqs_after)

    cat(sprintf("  BEFORE: n=%d | median=%d bp | range=[%d, %d]\n",
                length(len_before), median(len_before), min(len_before), max(len_before)))
    cat(sprintf("  AFTER:  n=%d | median=%d bp | range=[%d, %d]\n",
                length(len_after),  median(len_after),  min(len_after),  max(len_after)))

    # Coefficient of variation (CV = SD/mean) should decrease post-trim.
    # SeCAT trims to a consensus region, so length variance should shrink.
    cv_before <- sd(len_before) / mean(len_before)
    cv_after  <- sd(len_after)  / mean(len_after)
    cv_pass   <- cv_after <= cv_before
    cat(sprintf("  Length CV: %.4f → %.4f %s\n",
                cv_before, cv_after, if (cv_pass) "[PASS]" else "[NOTE: CV increased]"))

    add_result("ALL", "0A", "Sequence QC", "Length CV (after <= before)",
               cv_after, sprintf("<= %.4f (before)", cv_before), cv_pass)

    # Length distribution plot
    len_df <- bind_rows(
      data.frame(Length = len_before, Dataset = "Before"),
      data.frame(Length = len_after,  Dataset = "After")
    )
    p_len <- ggplot(len_df, aes(x = Length, fill = Dataset)) +
      geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
      scale_fill_manual(values = c("Before" = "steelblue", "After" = "coral")) +
      labs(title = "Sequence Length Distribution Before/After Trimming",
           x = "Sequence length (bp)", y = "Count") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(PLOTS_DIR, "00_seq_length_distribution.pdf"),
           p_len, width = 10, height = 5)
    cat("  ✓ Length distribution plot saved\n")

    # --- 0B: GC content preservation ---
    # GC content is a fundamental sequence property. Trimming should remove
    # flanking regions symmetrically, preserving the overall GC distribution.
    # A shift > 2 percentage points would suggest compositional bias.
    cat("\nTier 0B: GC content preservation\n")
    gc_before <- Biostrings::letterFrequency(seqs_before, letters = "GC",
                                             as.prob = TRUE)[, 1]
    gc_after  <- Biostrings::letterFrequency(seqs_after,  letters = "GC",
                                             as.prob = TRUE)[, 1]

    # Median GC should be stable; systematic shift would indicate bias
    gc_shift <- abs(median(gc_after) - median(gc_before))
    gc_pass  <- gc_shift < 0.02   # <2 percentage points = negligible
    cat(sprintf("  Median GC: %.3f → %.3f | shift = %.4f %s\n",
                median(gc_before), median(gc_after), gc_shift,
                if (gc_pass) "[PASS]" else "[NOTE: >2pp shift]"))

    add_result("ALL", "0B", "Sequence QC", "Median GC shift (absolute)",
               gc_shift, "< 0.02 (negligible bias)", gc_pass)

    gc_df <- bind_rows(
      data.frame(GC = gc_before, Dataset = "Before"),
      data.frame(GC = gc_after,  Dataset = "After")
    )
    p_gc <- ggplot(gc_df, aes(x = GC, fill = Dataset)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("Before" = "steelblue", "After" = "coral")) +
      labs(title = "GC Content Distribution Before/After Trimming",
           x = "GC proportion", y = "Density") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(PLOTS_DIR, "00b_gc_content.pdf"),
           p_gc, width = 8, height = 5)
    cat("  ✓ GC content plot saved\n")

    seq_qc_results$len_before  <- len_before
    seq_qc_results$len_after   <- len_after
    seq_qc_results$gc_before   <- gc_before
    seq_qc_results$gc_after    <- gc_after

  } else {
    cat("  ⚠ sequences.fasta not found in one or both directories — skipping Tier 0A/0B\n")
    cat("    (Place sequences.fasta in unaligned_cleaned/ and aligned_trimmed/ to enable)\n\n")
  }

  # --- 0C: Consensus coordinate verification ---
  # Cross-references the consensus region boundaries (from the alignment stage)
  # against the actual trimmed sequence lengths reported in trim_summary.csv.
  # Verifies that trimming hit the intended target region.
  cat("\nTier 0C: Consensus coordinate verification\n")
  if (file.exists(CONSENSUS_INFO) && file.exists(TRIM_SUMMARY)) {
    consensus_info <- read_csv(CONSENSUS_INFO, show_col_types = FALSE)
    trim_summary   <- read_csv(TRIM_SUMMARY,   show_col_types = FALSE)

    cons_start <- consensus_info$ConsensusStart[1]
    cons_end   <- consensus_info$ConsensusEnd[1]
    cons_len   <- cons_end - cons_start + 1

    cat(sprintf("  Consensus region: %d – %d (%d bp)\n", cons_start, cons_end, cons_len))
    cat(sprintf("  Studies processed: %d\n", nrow(trim_summary)))

    success_rate <- mean(trim_summary$status == "SUCCESS", na.rm = TRUE)
    success_pass <- success_rate >= 0.90
    cat(sprintf("  Trim success rate: %.1f%% %s\n",
                100 * success_rate,
                if (success_pass) "[PASS]" else "[NOTE: <90% success]"))

    if ("degapped_length_min" %in% colnames(trim_summary)) {
      ok_studies <- trim_summary %>% filter(status == "SUCCESS")
      # Trimmed sequences should be within ±10 bp of consensus length
      length_ok  <- ok_studies %>%
        mutate(ok = abs(degapped_length_min - cons_len) <= 10 |
                    abs(degapped_length_max - cons_len) <= 10)
      coord_pass <- mean(length_ok$ok) >= 0.80
      cat(sprintf("  Studies within ±10 bp of consensus length: %.1f%% %s\n",
                  100 * mean(length_ok$ok),
                  if (coord_pass) "[PASS]" else "[NOTE: coordinate mismatch]"))
      add_result("ALL", "0C", "Coordinate QC", "Studies within ±10 bp of target",
                 mean(length_ok$ok), ">= 0.80", coord_pass)
    }

    add_result("ALL", "0C", "Coordinate QC", "Trim success rate",
               success_rate, ">= 0.90", success_pass)

    write_csv(trim_summary %>%
                mutate(consensus_start = cons_start,
                       consensus_end   = cons_end,
                       consensus_len   = cons_len),
              file.path(OUTPUT_DIR, "tier0_coordinate_verification.csv"))
    cat("  ✓ Coordinate verification table saved\n")

  } else {
    cat("  ⚠ consensusregioninfo.csv or trim_summary.csv not found — skipping Tier 0C\n")
    cat("    (Run SeCAT pipeline through Stage 12 first)\n")
  }

  cat("\n✓ Tier 0 complete\n\n")

  # --- Load count data ---
  # Load pre- and post-trimming OTU tables, taxonomy, and metadata. These are
  # used for all subsequent tiers. Samples are aligned (intersected) to ensure
  # paired comparison; features are intersected with their taxonomy tables.
  cat("================================================================================\n")
  cat("LOADING COUNT DATA\n")
  cat("================================================================================\n\n")

  cat("Loading BEFORE data:\n")
  otu_before <- load_otu_table(PRE_OTU)
  tax_before <- load_taxonomy(PRE_TAX)
  meta_before <- load_metadata(PRE_META)

  cat("\nLoading AFTER data:\n")
  otu_after <- load_otu_table(POST_OTU)
  tax_after <- load_taxonomy(POST_TAX)
  meta_after <- load_metadata(POST_META)

  # Align samples -- only samples present in BOTH datasets can be compared.
  # This is essential for paired statistical tests (ICC, Bland-Altman, etc.).
  common_samples <- intersect(colnames(otu_before), colnames(otu_after))
  if (length(common_samples) == 0) stop("FATAL: No common samples between datasets.")
  cat(sprintf("\nCommon samples: %d\n", length(common_samples)))

  otu_before  <- otu_before[,  common_samples, drop = FALSE]
  otu_after   <- otu_after[,   common_samples, drop = FALSE]
  meta_before <- meta_before[common_samples,  , drop = FALSE]
  meta_after  <- meta_after[common_samples,   , drop = FALSE]

  # Align features to their respective taxonomy tables
  cf_before <- intersect(rownames(otu_before), rownames(tax_before))
  cf_after  <- intersect(rownames(otu_after),  rownames(tax_after))
  otu_before <- otu_before[cf_before, , drop = FALSE]
  tax_before <- tax_before[cf_before, , drop = FALSE]
  otu_after  <- otu_after[cf_after,   , drop = FALSE]
  tax_after  <- tax_after[cf_after,   , drop = FALSE]

  cat(sprintf("Features BEFORE: %d | AFTER: %d\n\n", nrow(otu_before), nrow(otu_after)))

  # Build ASV-level microtable objects (microeco package).
  # These serve as the base for all taxonomic levels -- aggregate_to_level()
  # collapses them to Genus or Family as needed in the validation loop.
  dataset_before_asv <- microtable$new(
    otu_table    = otu_before,
    tax_table    = tax_before,
    sample_table = meta_before
  )
  dataset_after_asv <- microtable$new(
    otu_table    = otu_after,
    tax_table    = tax_after,
    sample_table = meta_after
  )
  dataset_before_asv$tidy_dataset()
  dataset_after_asv$tidy_dataset()

  cat(sprintf("BEFORE (ASV): %d features × %d samples\n",
              nrow(dataset_before_asv$otu_table), ncol(dataset_before_asv$otu_table)))
  cat(sprintf("AFTER  (ASV): %d features × %d samples\n",
              nrow(dataset_after_asv$otu_table), ncol(dataset_after_asv$otu_table)))

  save_checkpoint("data", list(
    dataset_before_asv = dataset_before_asv,
    dataset_after_asv  = dataset_after_asv,
    meta_before        = meta_before,
    meta_after         = meta_after,
    otu_before         = otu_before,
    otu_after          = otu_after,
    tax_before         = tax_before,
    tax_after          = tax_after,
    summary_results    = summary_results,
    seq_qc_results     = seq_qc_results
  ))

} else {
  # Resume from data checkpoint
  ckpt <- load_checkpoint("data")
  if (is.null(ckpt)) stop("Cannot find data checkpoint. Run from START_TIER=0.")
  list2env(ckpt, envir = .GlobalEnv)
  cat(sprintf("Resumed: %d samples | BEFORE %d ASVs | AFTER %d ASVs\n\n",
              length(common_samples <- colnames(dataset_before_asv$otu_table)),
              nrow(dataset_before_asv$otu_table),
              nrow(dataset_after_asv$otu_table)))
}

# ==============================================================================
# MULTI-LEVEL VALIDATION LOOP
# ==============================================================================
# Iterates over each requested taxonomic resolution (default: ASV, Genus,
# Family). At each level, the full tier suite (1A through 4) is executed.
# This multi-resolution approach is important because:
#   - ASV level: tests whether individual sequence variants are preserved
#   - Genus level: tests whether biologically meaningful groups are stable
#   - Family level: tests robustness at coarser resolution (more conservative)
# If trimming preserves ecological signal at all three levels, it provides
# strong evidence that the harmonisation did not introduce systematic bias.
# ==============================================================================

for (taxa_level in TAXONOMIC_LEVELS) {

  cat(sprintf("\n################################################################################\n"))
  cat(sprintf("### TAXONOMIC LEVEL: %s ###\n", toupper(taxa_level)))
  cat(sprintf("################################################################################\n\n"))

  level_out  <- file.path(OUTPUT_DIR, tolower(taxa_level))
  level_plot <- file.path(PLOTS_DIR,  tolower(taxa_level))
  dir.create(level_out,  showWarnings = FALSE, recursive = TRUE)
  dir.create(level_plot, showWarnings = FALSE, recursive = TRUE)

  # Aggregate to this level
  cat("Preparing datasets...\n")
  dataset_before <- aggregate_to_level(dataset_before_asv, taxa_level)
  dataset_after  <- aggregate_to_level(dataset_after_asv,  taxa_level)
  dataset_before <- make_taxa_unique(dataset_before, taxa_level)
  dataset_after  <- make_taxa_unique(dataset_after,  taxa_level)

  # Compute alpha diversity indices (Observed, Chao1, ACE, Shannon, Simpson,
  # Pielou) for both datasets. PD (phylogenetic diversity) is disabled as it
  # requires a phylogenetic tree not available in the harmonised output.
  dataset_before$cal_alphadiv(PD = FALSE)
  dataset_after$cal_alphadiv(PD = FALSE)
  cat("  ✓ Alpha diversity calculated\n")

  # Compute beta diversity (Bray-Curtis dissimilarity). UniFrac is disabled
  # because it requires a phylogenetic tree. Bray-Curtis is the standard
  # distance metric for amplicon studies -- it is abundance-weighted and
  # bounded [0,1], making it suitable for community composition comparison.
  dataset_before$cal_betadiv(unifrac = FALSE)
  dataset_after$cal_betadiv(unifrac = FALSE)
  cat("  ✓ Beta diversity calculated\n\n")

  # ============================================================================
  # TIER 1A: ALPHA DIVERSITY PRESERVATION
  # ============================================================================
  # Tests whether per-sample diversity metrics are preserved after trimming.
  # Uses two complementary approaches:
  #
  #   ICC (Intraclass Correlation Coefficient, two-way agreement model):
  #     Measures absolute agreement between before/after values. ICC > 0.90
  #     indicates "excellent" agreement (Koo & Li, 2016). Unlike Pearson r,
  #     ICC penalises systematic offsets, not just rank changes.
  #
  #   Spearman rho: Rank-based correlation as a complementary check. Less
  #     sensitive to absolute scale shifts, but robust to outliers.
  #
  #   Bland-Altman plots: Visual assessment of agreement. Points should
  #     scatter symmetrically around zero bias with narrow limits of agreement.
  #
  # Metrics tested: Observed richness, Chao1 (bias-corrected richness),
  #   ACE (abundance-based coverage estimator), Shannon H' (entropy-based
  #   diversity), Simpson 1-D (dominance-weighted), Pielou J' (evenness).
  # ============================================================================

  cat(sprintf("--- TIER 1A: Alpha Diversity (%s) ---\n\n", taxa_level))

  alpha_b <- dataset_before$alpha_diversity
  alpha_a <- dataset_after$alpha_diversity

  common_alpha <- intersect(rownames(alpha_b), rownames(alpha_a))
  alpha_b <- alpha_b[common_alpha, , drop = FALSE]
  alpha_a <- alpha_a[common_alpha, , drop = FALSE]

  metrics <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "Pielou")
  metrics <- metrics[metrics %in% colnames(alpha_b) & metrics %in% colnames(alpha_a)]

  icc_results <- list()
  cor_results <- list()

  cat("ICC (two-way agreement) and Spearman rho:\n")
  for (m in metrics) {
    bv <- alpha_b[[m]]; av <- alpha_a[[m]]
    idx <- !is.na(bv) & !is.na(av)
    bv  <- bv[idx]; av <- av[idx]
    if (length(bv) < 3) next

    # Two-way random effects model with absolute agreement -- appropriate when
    # both raters (before/after) are the only raters of interest and we care
    # about agreement in actual values, not just consistency.
    icc_res <- irr::icc(data.frame(Before = bv, After = av),
                        model = "twoway", type = "agreement")
    icc_val <- icc_res$value
    icc_results[[m]] <- icc_val

    suppressWarnings({
      cor_res <- cor.test(bv, av, method = "spearman")
    })
    cor_results[[m]] <- cor_res$estimate

    pass <- icc_val > 0.90
    cat(sprintf("  %-12s ICC=%.4f | rho=%.4f (p=%.3e) %s\n",
                m, icc_val, cor_res$estimate, cor_res$p.value,
                if (pass) "[PASS]" else "[MODERATE]"))

    add_result(taxa_level, "1A", "Alpha Diversity",
               paste(m, "ICC"), icc_val, "> 0.90", pass)
  }

  # Scatter plots: before vs after for each diversity metric. The red dashed
  # 1:1 line shows perfect agreement; points above/below indicate systematic
  # increase/decrease. The blue regression line reveals any proportional bias.
  scatter_list <- list()
  for (m in metrics) {
    bv <- alpha_b[[m]]; av <- alpha_a[[m]]
    idx <- !is.na(bv) & !is.na(av)
    if (sum(idx) < 3) next
    df <- data.frame(Before = bv[idx], After = av[idx])
    p <- ggplot(df, aes(Before, After)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "navy") +
      annotate("text", x = -Inf, y = Inf,
               label = sprintf("ICC: %.3f\nrho: %.3f",
                               icc_results[[m]], cor_results[[m]]),
               hjust = -0.1, vjust = 1.2, size = 3.8) +
      labs(title = paste(m, "preservation"),
           x = paste(m, "(before)"), y = paste(m, "(after)")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    scatter_list[[m]] <- p
  }
  if (length(scatter_list) > 0) {
    ggsave(file.path(level_plot, sprintf("01_%s_alpha_scatter.pdf", tolower(taxa_level))),
           wrap_plots(scatter_list, ncol = 2), width = 12, height = 10)
    cat("  ✓ Alpha scatter plots saved\n")
  }

  # Bland-Altman plots
  ba_list <- list()
  for (m in metrics) {
    bv <- alpha_b[[m]]; av <- alpha_a[[m]]
    idx <- !is.na(bv) & !is.na(av)
    if (sum(idx) < 3) next
    ba_list[[m]] <- bland_altman_plot(bv[idx], av[idx],
                                     title = paste(m, "agreement"),
                                     xlab  = paste("Mean", m))
  }
  if (length(ba_list) > 0) {
    ggsave(file.path(level_plot, sprintf("02_%s_alpha_ba.pdf", tolower(taxa_level))),
           wrap_plots(ba_list, ncol = 2), width = 12, height = 10)
    cat("  ✓ Bland-Altman plots saved\n")
  }

  # ============================================================================
  # TIER 1B: BETA DIVERSITY PRESERVATION
  # ============================================================================
  # Tests whether the overall community structure (inter-sample relationships)
  # is preserved after trimming, using three complementary approaches:
  #
  #   Mantel test: Pearson correlation between the two Bray-Curtis distance
  #     matrices (before vs after). Tests whether samples that were similar
  #     before trimming remain similar after. r > 0.85 = strong preservation.
  #     Uses 999 permutations for significance testing.
  #
  #   Procrustes analysis: Rotates, reflects, and scales the post-trimming
  #     PCoA ordination to best fit the pre-trimming ordination. The residual
  #     sum of squares (m2) quantifies misfit. m2 < 0.10 = excellent fit.
  #     PROTEST adds a permutation test for the correlation between the two
  #     configurations.
  #
  #   Delta-R2 PERMANOVA: Compares the proportion of variance explained by
  #     StudyID before and after trimming. If trimming inflates the study
  #     effect (delta-R2 > 0.05), it may be introducing batch artefacts.
  #     Ideally delta-R2 <= 0 (study effect reduced or unchanged).
  # ============================================================================

  cat(sprintf("\n--- TIER 1B: Beta Diversity (%s) ---\n\n", taxa_level))

  dist_b <- as.matrix(dataset_before$beta_diversity$bray)
  dist_a <- as.matrix(dataset_after$beta_diversity$bray)

  common_dist <- intersect(rownames(dist_b), rownames(dist_a))
  db <- dist_b[common_dist, common_dist]
  da <- dist_a[common_dist, common_dist]

  # Mantel test: correlates the lower triangles of the two distance matrices.
  # Pearson method chosen over Spearman because Bray-Curtis distances are
  # bounded and approximately normally distributed for most community datasets.
  # 999 permutations: standard for ecological studies (Anderson, 2001).
  set.seed(42)
  mantel_res <- vegan::mantel(as.dist(db), as.dist(da),
                              method = "pearson", permutations = 999)
  mantel_pass <- mantel_res$statistic > 0.85
  cat(sprintf("  Mantel r = %.4f (p = %.4f) %s\n",
              mantel_res$statistic, mantel_res$signif,
              if (mantel_pass) "[PASS]" else "[MODERATE]"))
  add_result(taxa_level, "1B", "Beta Diversity", "Mantel r",
             mantel_res$statistic, "> 0.85", mantel_pass)

  # Procrustes / PROTEST analysis
  # PCoA (Principal Coordinates Analysis, aka classical MDS) embeds samples
  # in Euclidean space preserving Bray-Curtis distances. k=3 axes capture
  # the major variance dimensions. Procrustes then finds the optimal rotation
  # to superimpose the "after" ordination onto the "before" ordination.
  # symmetric=TRUE uses the symmetric Procrustes statistic (Peres-Neto & Jackson, 2001).
  pcoa_b <- cmdscale(as.dist(db), k = 3, eig = TRUE)
  pcoa_a <- cmdscale(as.dist(da), k = 3, eig = TRUE)
  proc   <- procrustes(pcoa_b$points, pcoa_a$points, symmetric = TRUE)
  set.seed(42)
  prot   <- protest(pcoa_b$points, pcoa_a$points, permutations = 999)

  proc_pass <- proc$ss  < 0.10
  prot_pass <- prot$t0  > 0.95
  cat(sprintf("  Procrustes m² = %.4f %s\n",
              proc$ss, if (proc_pass) "[PASS]" else "[MODERATE]"))
  cat(sprintf("  PROTEST r = %.4f (p = %.4f) %s\n",
              prot$t0, prot$signif, if (prot_pass) "[PASS]" else "[MODERATE]"))
  add_result(taxa_level, "1B", "Beta Diversity", "Procrustes m²",
             proc$ss, "< 0.10", proc_pass)
  add_result(taxa_level, "1B", "Beta Diversity", "PROTEST correlation",
             prot$t0, "> 0.95", prot_pass)

  # PERMANOVA (adonis2) -- delta-R2 framing
  # PERMANOVA partitions variance in the distance matrix by a grouping factor
  # (StudyID). The key metric is the CHANGE in R2, not the R2 itself.
  # Assumptions: PERMANOVA is sensitive to differences in dispersion between
  # groups (Anderson, 2001); betadisper in Tier 1C checks this assumption.
  # A negative delta-R2 means trimming reduced the study batch effect (good).
  has_study_id <- "StudyID" %in% colnames(meta_before)
  if (has_study_id) {
    tryCatch({
      meta_dist <- meta_before[common_dist, , drop = FALSE]
      perm_b <- vegan::adonis2(as.dist(db) ~ StudyID, data = meta_dist,
                               permutations = 999)
      perm_a <- vegan::adonis2(as.dist(da) ~ StudyID, data = meta_dist,
                               permutations = 999)
      r2_b   <- perm_b$R2[1]
      r2_a   <- perm_a$R2[1]
      dr2    <- r2_a - r2_b  # negative = batch effect reduced
      dr2_pass <- dr2 <= 0.05   # trimming should not inflate study effect

      cat(sprintf("  PERMANOVA R²: %.3f → %.3f | ΔR² = %+.3f %s\n",
                  r2_b, r2_a, dr2,
                  if (dr2_pass) "[PASS]" else "[NOTE: study effect increased]"))
      add_result(taxa_level, "1B", "PERMANOVA", "ΔR² (after − before)",
                 dr2, "<= 0.05", dr2_pass)

      write_csv(bind_rows(
        as_tibble(perm_b, rownames = "Term") %>% mutate(Dataset = "Before"),
        as_tibble(perm_a, rownames = "Term") %>% mutate(Dataset = "After")
      ), file.path(level_out, "permanova_results.csv"))
      cat("  ✓ PERMANOVA results saved\n")
    }, error = function(e) cat("  ⚠ PERMANOVA failed:", conditionMessage(e), "\n"))
  }

  # PCoA ordination plots (before vs after)
  # Uses microeco's trans_beta class for consistent plot styling. Ellipses
  # (convex hulls) are drawn per StudyID when available, allowing visual
  # assessment of whether study clusters maintain their separation pattern.
  tryCatch({
    beta_b <- trans_beta$new(
      dataset = dataset_before,
      measure = "bray",
      group   = if (has_study_id) "StudyID" else NULL
    )
    beta_a <- trans_beta$new(
      dataset = dataset_after,
      measure = "bray",
      group   = if (has_study_id) "StudyID" else NULL
    )
    beta_b$cal_ordination(method = "PCoA")
    beta_a$cal_ordination(method = "PCoA")

    plot_args <- list(plot_type = if (has_study_id) c("point", "ellipse") else "point")
    if (has_study_id) {
      plot_args$plot_color <- "StudyID"
      plot_args$plot_shape <- "StudyID"
      plot_args$ellipse_chull_fill <- FALSE
    }

    p_b <- do.call(beta_b$plot_ordination, plot_args) +
      labs(title = "PCoA before trimming") + theme_bw(base_size = 11)
    p_a <- do.call(beta_a$plot_ordination, plot_args) +
      labs(title = "PCoA after trimming") + theme_bw(base_size = 11)

    ggsave(file.path(level_plot, sprintf("04_%s_pcoa_before.pdf", tolower(taxa_level))),
           p_b, width = 10, height = 8)
    ggsave(file.path(level_plot, sprintf("05_%s_pcoa_after.pdf", tolower(taxa_level))),
           p_a, width = 10, height = 8)
    cat("  ✓ PCoA plots saved\n")
  }, error = function(e) cat("  ⚠ PCoA plots failed:", conditionMessage(e), "\n"))

  # Procrustes overlay plot -- shows each sample as a pair of points (before
  # and after) connected by an arrow. Short arrows = good preservation.
  # Colour-coded by StudyID to reveal whether any study is disproportionately
  # affected by trimming.
  tryCatch({
    proc_df <- bind_rows(
      data.frame(PC1 = pcoa_b$points[, 1], PC2 = pcoa_b$points[, 2],
                 SampleID = rownames(pcoa_b$points), Dataset = "Before"),
      data.frame(PC1 = proc$Yrot[, 1],     PC2 = proc$Yrot[, 2],
                 SampleID = rownames(proc$Yrot),     Dataset = "After")
    )
    arrow_df <- data.frame(
      x = pcoa_b$points[, 1], y = pcoa_b$points[, 2],
      xend = proc$Yrot[, 1],  yend = proc$Yrot[, 2]
    )

    if (has_study_id) {
      meta_sub <- meta_before %>%
        rownames_to_column("SampleID") %>%
        filter(SampleID %in% rownames(pcoa_b$points)) %>%
        select(SampleID, StudyID)
      proc_df <- left_join(proc_df, meta_sub, by = "SampleID")

      p_proc <- ggplot() +
        geom_segment(data = arrow_df,
                     aes(x=x, y=y, xend=xend, yend=yend),
                     arrow = arrow(length = unit(0.012, "npc")),
                     alpha = 0.25, color = "grey40", linewidth = 0.3) +
        geom_point(data = proc_df,
                   aes(PC1, PC2, color = StudyID, shape = Dataset),
                   size = 3, alpha = 0.8) +
        scale_shape_manual(values = c("Before" = 16, "After" = 17)) +
        scale_color_brewer(palette = "Set2", name = "Study") +
        labs(title = "Procrustes overlay",
             subtitle = sprintf("m² = %.4f | r = %.4f (p = %.4f)",
                                proc$ss, prot$t0, prot$signif),
             x = "PC1", y = "PC2") +
        theme_bw(base_size = 11) +
        theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 9))
    } else {
      p_proc <- ggplot() +
        geom_segment(data = arrow_df,
                     aes(x=x, y=y, xend=xend, yend=yend),
                     arrow = arrow(length = unit(0.012, "npc")),
                     alpha = 0.25, color = "grey40", linewidth = 0.3) +
        geom_point(data = proc_df,
                   aes(PC1, PC2, color = Dataset, shape = Dataset),
                   size = 3, alpha = 0.8) +
        scale_shape_manual(values = c("Before" = 16, "After" = 17)) +
        labs(title = "Procrustes overlay",
             subtitle = sprintf("m² = %.4f | r = %.4f",
                                proc$ss, prot$t0),
             x = "PC1", y = "PC2") +
        theme_bw(base_size = 11) +
        theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 9))
    }

    ggsave(file.path(level_plot, sprintf("06_%s_procrustes.pdf", tolower(taxa_level))),
           p_proc, width = 10, height = 8)
    cat("  ✓ Procrustes overlay saved\n")
  }, error = function(e) cat("  ⚠ Procrustes plot failed:", conditionMessage(e), "\n"))

  # ============================================================================
  # TIER 1C: WITHIN-GROUP DISPERSION (betadisper / Levene analogue)
  # ============================================================================
  # betadisper (Anderson, 2006) is the multivariate analogue of Levene's test.
  # It computes each sample's distance to its group centroid in ordination
  # space, then tests whether these distances differ between groups.
  #
  # Why this matters for validation:
  #   - PERMANOVA is sensitive to heterogeneous dispersion (confounding).
  #   - If trimming systematically shrinks or inflates within-study variance,
  #     it would change the PERMANOVA interpretation.
  #   - A relative change < 10% in mean distance-to-centroid is considered
  #     negligible for this validation.
  #
  # The permutation test (permutest) assesses significance of dispersion
  # differences among studies, in both the before and after datasets.
  # ============================================================================

  if (has_study_id) {
    cat(sprintf("\n--- TIER 1C: Within-group Dispersion (%s) ---\n\n", taxa_level))
    tryCatch({
      meta_dist <- meta_before[common_dist, , drop = FALSE]
      grp <- as.factor(meta_dist$StudyID)

      disp_b <- vegan::betadisper(as.dist(db), grp)
      disp_a <- vegan::betadisper(as.dist(da), grp)

      # Average distance to centroid per group
      avg_disp_b <- mean(disp_b$distances, na.rm = TRUE)
      avg_disp_a <- mean(disp_a$distances, na.rm = TRUE)
      disp_change <- abs(avg_disp_a - avg_disp_b) / avg_disp_b

      # Permutation test for difference in dispersion
      set.seed(42)
      lev_b <- vegan::permutest(disp_b, permutations = 999)
      lev_a <- vegan::permutest(disp_a, permutations = 999)

      disp_pass <- disp_change < 0.10  # <10% relative change in avg dispersion
      cat(sprintf("  Mean distance-to-centroid: %.4f → %.4f | relative change = %.1f%% %s\n",
                  avg_disp_b, avg_disp_a, 100 * disp_change,
                  if (disp_pass) "[PASS]" else "[NOTE: >10% change]"))
      add_result(taxa_level, "1C", "Dispersion", "Relative change in mean dispersion",
                 disp_change, "< 0.10 (< 10% change)", disp_pass)

      # Dispersion plot
      disp_df <- bind_rows(
        data.frame(Distance = disp_b$distances, Group = grp, Dataset = "Before"),
        data.frame(Distance = disp_a$distances, Group = grp, Dataset = "After")
      )
      p_disp <- ggplot(disp_df, aes(x = Group, y = Distance, fill = Dataset)) +
        geom_boxplot(alpha = 0.6, outlier.size = 1.5) +
        scale_fill_manual(values = c("Before" = "steelblue", "After" = "coral")) +
        labs(title = "Within-group dispersion (distance to centroid)",
             x = "Study", y = "Distance to centroid") +
        theme_bw(base_size = 11) +
        theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
              axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(level_plot, sprintf("07_%s_dispersion.pdf", tolower(taxa_level))),
             p_disp, width = 10, height = 6)
      cat("  ✓ Dispersion plot saved\n")

    }, error = function(e) cat("  ⚠ Dispersion test failed:", conditionMessage(e), "\n"))
  } else {
    cat("  ⚠ StudyID not found — Tier 1C skipped\n")
  }

  # ============================================================================
  # TIER 1D: RANK ABUNDANCE CURVE PRESERVATION
  # ============================================================================
  # Rank abundance (Whittaker) curves describe the dominance structure of a
  # community: from the most abundant taxon to the rarest. This tier tests
  # whether the shape of this curve is preserved after trimming.
  #
  # For each sample, relative abundances are sorted in decreasing order for
  # both before and after datasets, then Spearman rho is computed between the
  # two sorted vectors. High rho (> 0.90) = the dominant-to-rare gradient
  # is preserved. This is important because ecological conclusions about
  # community evenness and dominance depend on this structure.
  # ============================================================================

  cat(sprintf("\n--- TIER 1D: Rank Abundance Preservation (%s) ---\n\n", taxa_level))
  tryCatch({
    otu_b_mat <- as.matrix(dataset_before$otu_table)
    otu_a_mat <- as.matrix(dataset_after$otu_table)

    # Relative abundances per sample (column-wise normalisation to proportions).
    # This removes library size effects, focusing on compositional structure.
    rel_b <- sweep(otu_b_mat, 2, colSums(otu_b_mat), "/")
    rel_a <- sweep(otu_a_mat, 2, colSums(otu_a_mat), "/")

    common_ra <- intersect(colnames(rel_b), colnames(rel_a))

    ra_cors <- sapply(common_ra, function(s) {
      rv_b <- sort(rel_b[, s], decreasing = TRUE)
      rv_a <- sort(rel_a[, s], decreasing = TRUE)
      # Compare top-N shared length
      n <- min(length(rv_b), length(rv_a))
      if (n < 3) return(NA_real_)
      suppressWarnings(cor(rv_b[1:n], rv_a[1:n], method = "spearman"))
    })

    mean_ra_cor  <- mean(ra_cors, na.rm = TRUE)
    prop_ra_pass <- mean(ra_cors >= 0.90, na.rm = TRUE)
    ra_pass      <- mean_ra_cor > 0.90

    cat(sprintf("  Mean Spearman rho (rank abundance): %.4f %s\n",
                mean_ra_cor, if (ra_pass) "[PASS]" else "[MODERATE]"))
    cat(sprintf("  Proportion samples with rho >= 0.90: %.1f%%\n",
                100 * prop_ra_pass))

    add_result(taxa_level, "1D", "Rank Abundance",
               "Mean Spearman rho", mean_ra_cor, "> 0.90", ra_pass)
    add_result(taxa_level, "1D", "Rank Abundance",
               "Proportion samples rho >= 0.90", prop_ra_pass, "> 0.80",
               prop_ra_pass > 0.80)

    # Distribution of per-sample rho
    p_ra <- ggplot(data.frame(rho = ra_cors), aes(x = rho)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
      geom_vline(xintercept = mean_ra_cor, color = "red",
                 linetype = "dashed", linewidth = 0.9) +
      annotate("text", x = mean_ra_cor, y = Inf,
               label = sprintf("Mean rho = %.3f", mean_ra_cor),
               hjust = -0.1, vjust = 1.5, color = "red", size = 3.5) +
      labs(title = "Per-sample rank abundance curve preservation",
           x = "Spearman rho (before vs after)", y = "Samples") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(level_plot, sprintf("08_%s_rank_abundance.pdf", tolower(taxa_level))),
           p_ra, width = 8, height = 5)
    cat("  ✓ Rank abundance plot saved\n")

  }, error = function(e) cat("  ⚠ Rank abundance test failed:", conditionMessage(e), "\n"))

  # ============================================================================
  # TIER 2A: CROSS-STUDY FEATURE SHARING
  # ============================================================================
  # A core goal of SeCAT is to enable cross-study comparability. This tier
  # tests whether trimming improves (or at least maintains) the degree of
  # feature sharing between studies.
  #
  # At ASV level: Uses the MetaASV mapping to determine how many consensus
  #   ASVs (MetaASVs) are detected in >= 2 studies. The rate should increase
  #   or remain stable after trimming, because trimming to a common region
  #   allows identical sequences from different studies to be recognised as
  #   the same ASV.
  #
  # At Genus/Family level: Computes pairwise Jaccard similarity (intersection
  #   over union) between the taxa detected in each study pair. A stable or
  #   increasing mean Jaccard indicates that trimming did not cause
  #   differential taxon dropout across studies.

  cat(sprintf("\n--- TIER 2A: Cross-study Sharing (%s) ---\n\n", taxa_level))

  if (has_study_id) {
    studies <- unique(meta_before$StudyID)
    studies <- studies[!is.na(studies)]          # drop NA study IDs before matrix build
    n_studies <- length(studies)
    cat(sprintf("  Studies: %d (%s)\n", n_studies, paste(studies, collapse=", ")))

    if (taxa_level == "ASV" && file.exists(ASV_MAPPING)) {
      # -----------------------------------------------------------------------
      # ASV level: MetaASV sharing rate
      # Uses asv_mapping_final.tsv to bridge original ASV hashes to MetaASV
      # IDs. For each MetaASV, counts how many distinct studies contribute
      # detections. The sharing rate = proportion of MetaASVs in >= 2 studies.
      # -----------------------------------------------------------------------
      cat("  Using MetaASV mapping for ASV-level cross-study sharing\n")
      asv_map <- read_tsv(ASV_MAPPING, show_col_types = FALSE)

      # Attach study membership from metadata
      # Assumes mapping has columns: Original_ID (or Hash), Meta_ID
      # and OTU table columns are sample IDs
      otu_b_df <- dataset_before_asv$otu_table %>%
        rownames_to_column("Hash")
      meta_b_df <- meta_before %>% rownames_to_column("SampleID")

      # For each MetaASV, count how many studies it appears in (after mapping)
      hash_to_meta <- asv_map %>%
        mutate(Hash = str_extract(Original_ID, "[0-9a-f]{32}$")) %>%
        filter(!is.na(Hash)) %>%
        distinct(Hash, Meta_ID)

      study_presence <- otu_b_df %>%
        left_join(hash_to_meta, by = "Hash") %>%
        filter(!is.na(Meta_ID)) %>%
        pivot_longer(cols = -c(Hash, Meta_ID), names_to = "SampleID",
                     values_to = "Count") %>%
        filter(Count > 0) %>%
        left_join(meta_b_df %>% select(SampleID, StudyID), by = "SampleID") %>%
        group_by(Meta_ID) %>%
        summarise(n_studies = n_distinct(StudyID), .groups = "drop")

      shared_rate_before <- mean(study_presence$n_studies >= 2, na.rm = TRUE)

      # After: MetaASVs are the direct rownames of otu_after
      otu_a_df <- dataset_after_asv$otu_table %>%
        rownames_to_column("Meta_ID")
      meta_a_df <- meta_after %>% rownames_to_column("SampleID")

      study_presence_after <- otu_a_df %>%
        pivot_longer(cols = -Meta_ID, names_to = "SampleID", values_to = "Count") %>%
        filter(Count > 0) %>%
        left_join(meta_a_df %>% select(SampleID, StudyID), by = "SampleID") %>%
        group_by(Meta_ID) %>%
        summarise(n_studies = n_distinct(StudyID), .groups = "drop")

      shared_rate_after <- mean(study_presence_after$n_studies >= 2, na.rm = TRUE)

      share_pass <- shared_rate_after >= shared_rate_before ||
              abs(shared_rate_after - shared_rate_before) < 0.005
      cat(sprintf("  MetaASVs shared across ≥2 studies: %.1f%% → %.1f%% %s\n",
                  100 * shared_rate_before, 100 * shared_rate_after,
                  if (share_pass) "[PASS: improved/stable]" else "[NOTE: decreased]"))
      add_result(taxa_level, "2A", "Cross-study Sharing",
                 "MetaASV sharing rate (>= 2 studies)",
                 shared_rate_after, ">= before", share_pass)

    } else if (taxa_level %in% c("Genus", "Family")) {
      # -----------------------------------------------------------------------
      # Genus/Family level: pairwise Jaccard similarity
      # Jaccard index J(A,B) = |A intersect B| / |A union B|, where A and B
      # are the sets of taxa detected in two studies. Ranges from 0 (no
      # shared taxa) to 1 (identical taxa sets). Computed for all study pairs
      # both before and after trimming.
      # -----------------------------------------------------------------------
      # Jaccard index helper (set-based, presence/absence only)
      jaccard_idx <- function(a, b) {
        length(intersect(a, b)) / length(union(a, b))
      }

      taxa_by_study <- function(dataset, meta_df, level) {
        otu <- dataset$otu_table
        lapply(studies, function(st) {
          samps <- rownames(meta_df)[meta_df$StudyID == st]
          samps <- intersect(samps, colnames(otu))
          if (length(samps) == 0) return(character(0))
          rownames(otu)[rowSums(otu[, samps, drop = FALSE]) > 0]
        }) %>% setNames(studies)
      }

      tb_before <- taxa_by_study(dataset_before, meta_before, taxa_level)
      tb_after  <- taxa_by_study(dataset_after,  meta_after,  taxa_level)

      ji_mat_b <- ji_mat_a <- matrix(NA, n_studies, n_studies,
                                     dimnames = list(studies, studies))
      for (i in seq_len(n_studies)) {
        for (j in i:n_studies) {
          if (i == j) {
            ji_mat_b[i, j] <- ji_mat_a[i, j] <- 1
          } else {
            ji_mat_b[i, j] <- ji_mat_b[j, i] <-
              jaccard_idx(tb_before[[studies[i]]], tb_before[[studies[j]]])
            ji_mat_a[i, j] <- ji_mat_a[j, i] <-
              jaccard_idx(tb_after[[studies[i]]], tb_after[[studies[j]]])
          }
        }
      }

      mj_b <- mean(ji_mat_b[lower.tri(ji_mat_b)], na.rm = TRUE)
      mj_a <- mean(ji_mat_a[lower.tri(ji_mat_a)], na.rm = TRUE)
      ji_pass <- (mj_a >= mj_b) || (abs(mj_a - mj_b) / max(mj_b, 1e-9) < 0.02)
      change  <- 100 * (mj_a - mj_b) / max(mj_b, 1e-9)

      cat(sprintf("  Mean pairwise Jaccard: %.4f → %.4f (%+.1f%%) %s\n",
                  mj_b, mj_a, change,
                  if (ji_pass) "[PASS: improved/stable]" else "[NOTE: decreased]"))
      add_result(taxa_level, "2A", "Cross-study Sharing",
                 paste("Mean pairwise Jaccard", taxa_level),
                 mj_a, ">= before", ji_pass)

    write_csv(as.data.frame(ji_mat_b) %>% 
                `rownames<-`(make.names(rownames(ji_mat_b), unique=TRUE)) %>%
                rownames_to_column("Study"),
               file.path(level_out, "jaccard_before.csv"))
    write_csv(as.data.frame(ji_mat_a) %>% 
                `rownames<-`(make.names(rownames(ji_mat_a), unique=TRUE)) %>%
                rownames_to_column("Study"),
              file.path(level_out, "jaccard_after.csv"))
      cat("  ✓ Jaccard matrices saved\n")

    } else {
      cat("  ⚠ ASV mapping file not found — ASV-level cross-study sharing skipped\n")
      cat("    (Provide asv_mapping_final.tsv in BASE_DIR to enable)\n")
    }
  } else {
    cat("  ⚠ StudyID not in metadata — Tier 2A skipped\n")
  }

  # ============================================================================
  # TIER 3A: TAXONOMIC COMPOSITION BAR PLOTS
  # ============================================================================
cat(sprintf("--- TIER 3A: Taxonomic Composition (%s) ---\n\n", taxa_level))
display_rank <- if (taxa_level == "ASV") "Phylum" else taxa_level
tryCatch({
  tab_b <- trans_abund$new(dataset = dataset_before, taxrank = display_rank, ntaxa = 10,
                           groupmean = if (has_study_id) "StudyID" else NULL)
  tab_a <- trans_abund$new(dataset = dataset_after,  taxrank = display_rank, ntaxa = 10,
                           groupmean = if (has_study_id) "StudyID" else NULL)

  # Extract taxa names from the transformed data (data_abund, not abund_list)
  taxa_col <- colnames(tab_b$data_abund)[1]
  all_taxa <- union(unique(tab_b$data_abund[[taxa_col]]),
                    unique(tab_a$data_abund[[taxa_col]]))
  n_taxa   <- length(all_taxa)
  pal_cols <- if (n_taxa <= 12) {
    RColorBrewer::brewer.pal(max(3, n_taxa), "Paired")[seq_len(n_taxa)]
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_taxa)
  }
  shared_pal <- setNames(pal_cols, all_taxa)

  # Force consistent factor levels so stacking order matches between plots
  tab_b$data_abund[[taxa_col]] <- factor(tab_b$data_abund[[taxa_col]], levels = names(shared_pal))
  tab_a$data_abund[[taxa_col]] <- factor(tab_a$data_abund[[taxa_col]], levels = names(shared_pal))

  p_bar_b <- tab_b$plot_bar(legend_text_italic = FALSE) +
    scale_fill_manual(values = shared_pal, drop = FALSE) +
    labs(title = sprintf("%s composition before (%s level)", display_rank, taxa_level))
  p_bar_a <- tab_a$plot_bar(legend_text_italic = FALSE) +
    scale_fill_manual(values = shared_pal, drop = FALSE) +
    labs(title = sprintf("%s composition after (%s level)", display_rank, taxa_level))

  ggsave(file.path(level_plot, sprintf("09_%s_tax_before.pdf", tolower(taxa_level))),
         p_bar_b, width = 12, height = 7)
  ggsave(file.path(level_plot, sprintf("10_%s_tax_after.pdf",  tolower(taxa_level))),
         p_bar_a, width = 12, height = 7)
  cat("  ✓ Composition bar plots saved\n")
}, error = function(e) cat("  Composition plots failed:", conditionMessage(e), "\n"))

  # ============================================================================
  # TIER 3B: CO-OCCURRENCE NETWORK STABILITY
  # ============================================================================

  cat(sprintf("\n--- TIER 3B: Network Stability (%s) ---\n\n", taxa_level))

  tryCatch({
    cat("  Building co-occurrence networks (this may take several minutes)...\n")

    net_b <- trans_network$new(
      dataset     = dataset_before,
      cor_method  = "spearman",
      filter_thres = 0.0005,
      nThreads    = min(4, ncores)
    )
    net_b$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
    cat("  ✓ BEFORE network built\n")

    net_a <- trans_network$new(
      dataset     = dataset_after,
      cor_method  = "spearman",
      filter_thres = 0.0005,
      nThreads    = min(4, ncores)
    )
    net_a$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
    cat("  ✓ AFTER network built\n\n")

    g_b <- net_b$res_network
    g_a <- net_a$res_network

    verts_b <- igraph::vcount(g_b); verts_a <- igraph::vcount(g_a)
    edges_b <- igraph::ecount(g_b); edges_a <- igraph::ecount(g_a)
    vert_pct <- 100 * (verts_a - verts_b) / max(verts_b, 1)
    edge_pct <- 100 * (edges_a - edges_b) / max(edges_b, 1)

    cat(sprintf("  Vertices: %d → %d (%+.1f%%)\n", verts_b, verts_a, vert_pct))
    cat(sprintf("  Edges:    %d → %d (%+.1f%%)\n", edges_b, edges_a, edge_pct))

    add_result(taxa_level, "3B", "Network Stability",
               "Vertex count % change", vert_pct, "< 20%", abs(vert_pct) < 20)
    add_result(taxa_level, "3B", "Network Stability",
               "Edge count % change",   edge_pct, "< 20%", abs(edge_pct) < 20)

    # ------------------------------------------------------------------
    # Hub node preservation
    # Hub nodes = top 10% by degree in the BEFORE network.
    # At ASV level: BEFORE nodes are original hashes, AFTER nodes are
    # MetaASV IDs — names never match directly. Translate hub hashes
    # to MetaASV IDs via asv_mapping_final.tsv before comparing.
    # At Genus/Family level: names match directly.
    # ------------------------------------------------------------------
    deg_b   <- igraph::degree(g_b)
    hub_thr <- quantile(deg_b, 0.90)
    hubs_b  <- names(deg_b[deg_b >= hub_thr])
    nodes_a <- igraph::V(g_a)$name

    if (taxa_level == "ASV" && file.exists(ASV_MAPPING)) {
      hub_map <- read_tsv(ASV_MAPPING, show_col_types = FALSE) %>%
        mutate(Hash = str_extract(Original_ID, "[0-9a-f]{32}$")) %>%
        filter(!is.na(Hash)) %>%
        distinct(Hash, Meta_ID)
      # Translate hub hashes → MetaASV IDs; keep only those that mapped
      hubs_b_translated <- hub_map$Meta_ID[hub_map$Hash %in% hubs_b]
      n_mappable        <- length(hubs_b_translated)
      if (n_mappable == 0) {
        cat("  ⚠ Hub nodes: none of the hub hashes found in mapping file\n")
        hub_retained <- NA_real_
      } else {
        hub_retained <- mean(hubs_b_translated %in% nodes_a, na.rm = TRUE)
      }
      cat(sprintf("  Hub nodes retained (top 10%% degree, via MetaASV map): %d/%d (%.1f%%) %s\n",
                  sum(hubs_b_translated %in% nodes_a), n_mappable,
                  100 * ifelse(is.na(hub_retained), 0, hub_retained),
                  if (!is.na(hub_retained) && hub_retained >= 0.80) "[PASS]"
                  else "[NOTE: hub loss > 20%]"))
    } else {
      # Genus / Family: names match directly
      hub_retained <- mean(hubs_b %in% nodes_a, na.rm = TRUE)
      cat(sprintf("  Hub nodes retained (top 10%% degree): %d/%d (%.1f%%) %s\n",
                  sum(hubs_b %in% nodes_a), length(hubs_b),
                  100 * hub_retained,
                  if (hub_retained >= 0.80) "[PASS]" else "[NOTE: hub loss > 20%]"))
    }

    hub_pass <- !is.na(hub_retained) && hub_retained >= 0.80
    add_result(taxa_level, "3B", "Network Stability",
               "Hub node retention", ifelse(is.na(hub_retained), 0, hub_retained),
               ">= 0.80", hub_pass)

    # Network attribute table
    net_attrs <- tibble(
      Dataset      = c("Before", "After"),
      Taxa_level   = taxa_level,
      Vertices     = c(verts_b, verts_a),
      Edges        = c(edges_b, edges_a),
      Density      = c(igraph::edge_density(g_b), igraph::edge_density(g_a)),
      Mean_degree  = c(mean(igraph::degree(g_b)), mean(igraph::degree(g_a))),
      Transitivity = c(igraph::transitivity(g_b, type = "global"),
                       igraph::transitivity(g_a, type = "global")),
      Hub_retention = c(NA, hub_retained)
    )
    write_csv(net_attrs, file.path(level_out, "network_attributes.csv"))
    cat("  ✓ Network attributes saved\n")

    # meconetcomp analysis (if available)
    if (HAS_MECONETCOMP) {
      tryCatch({
        net_list <- list(Before = net_b, After = net_a)
        net_list <- cal_module(net_list, undirected_method = "cluster_fast_greedy")
        na_full  <- cal_network_attr(net_list)
        write_csv(na_full, file.path(level_out, "network_attributes_full.csv"))

        node_ov  <- node_comp(net_list, property = "name")
        node_vn  <- trans_venn$new(node_ov, ratio = "numratio")
        g_nv     <- node_vn$plot_venn(fill_color = FALSE)
        ggsave(file.path(level_plot, sprintf("14_%s_node_overlap.pdf", tolower(taxa_level))),
               g_nv, width = 7, height = 6)

        edge_ov  <- edge_comp(net_list)
        edge_vn  <- trans_venn$new(edge_ov, ratio = "numratio")
        g_ev     <- edge_vn$plot_venn(fill_color = FALSE)
        ggsave(file.path(level_plot, sprintf("15_%s_edge_overlap.pdf", tolower(taxa_level))),
               g_ev, width = 7, height = 6)

        rob <- robustness$new(net_list,
                              remove_strategy = c("edge_rand", "node_degree_high"),
                              remove_ratio    = seq(0, 0.99, 0.1),
                              measure         = "Eff", run = 10)
        ggsave(file.path(level_plot, sprintf("16_%s_robustness.pdf", tolower(taxa_level))),
               rob$plot(linewidth = 1), width = 12, height = 8)
        cat("  ✓ meconetcomp analyses saved\n")
      }, error = function(e) {
        cat("  ⚠ meconetcomp analysis failed:", conditionMessage(e), "\n")
      })
    }

# Export networks for Gephi
tryCatch({
  igraph::write_graph(g_b,
    file.path(level_out, sprintf("12_%s_network_before.graphml", tolower(taxa_level))),
    format = "graphml")
  igraph::write_graph(g_a,
    file.path(level_out, sprintf("13_%s_network_after.graphml", tolower(taxa_level))),
    format = "graphml")
  cat("  ✓ Network graphs exported for Gephi (GraphML format)\n")
}, error = function(e) cat("  ⚠ Network export failed:", conditionMessage(e), "\n"))

}, error = function(e) {
  cat(sprintf("  ⚠ Network analysis failed: %s\n", conditionMessage(e)))
})

  # ============================================================================
  # TIER 4: ABUNDANCE CONCORDANCE
  # ============================================================================
  # Replaces DESeq2. For each taxon present in both datasets, compute Spearman
  # correlation of abundance vectors across samples. Reports:
  #   - Mean per-taxon rho
  #   - Proportion of taxa with rho >= 0.90
  # At ASV level, uses MetaASV bridge; at Genus/Family, taxa names match directly.

  cat(sprintf("\n--- TIER 4: Abundance Concordance (%s) ---\n\n", taxa_level))

  tryCatch({
    otu_b_mat <- as.matrix(dataset_before$otu_table)
    otu_a_mat <- as.matrix(dataset_after$otu_table)

    # Relative abundances (per-sample normalisation avoids library size confound)
    rel_b <- sweep(otu_b_mat, 2, colSums(otu_b_mat), "/")
    rel_a <- sweep(otu_a_mat, 2, colSums(otu_a_mat), "/")

    if (taxa_level == "ASV" && file.exists(ASV_MAPPING)) {
      # Map BEFORE hashes to MetaASV IDs
      asv_map    <- read_tsv(ASV_MAPPING, show_col_types = FALSE)
      hash_map   <- asv_map %>%
        mutate(Hash = str_extract(Original_ID, "[0-9a-f]{32}$")) %>%
        filter(!is.na(Hash)) %>%
        distinct(Hash, Meta_ID)

      rel_b_df   <- as.data.frame(rel_b) %>% rownames_to_column("Hash")
      rel_b_meta <- left_join(rel_b_df, hash_map, by = "Hash") %>%
        filter(!is.na(Meta_ID)) %>%
        select(-Hash) %>%
        group_by(Meta_ID) %>%
        summarise(across(everything(), sum), .groups = "drop") %>%
        column_to_rownames("Meta_ID")

      common_taxa    <- intersect(rownames(rel_b_meta), rownames(rel_a))
      common_samples <- intersect(colnames(rel_b_meta), colnames(rel_a))

    } else {
      # Genus / Family: names match directly
      common_taxa    <- intersect(rownames(rel_b), rownames(rel_a))
      common_samples <- intersect(colnames(rel_b), colnames(rel_a))
      rel_b_meta     <- as.data.frame(rel_b)
    }

    cat(sprintf("  Aligned: %d taxa × %d samples\n",
                length(common_taxa), length(common_samples)))

    if (length(common_taxa) >= 5 && length(common_samples) >= 5) {
      rb_sub <- as.matrix(rel_b_meta[common_taxa, common_samples])
      ra_sub <- as.matrix(rel_a[common_taxa,    common_samples])

      # Per-taxon Spearman rho across samples
      taxon_cors <- sapply(common_taxa, function(tx) {
        suppressWarnings(
          cor(rb_sub[tx, ], ra_sub[tx, ], method = "spearman")
        )
      })

      mean_cor  <- mean(taxon_cors, na.rm = TRUE)
      prop_pass <- mean(taxon_cors >= 0.90, na.rm = TRUE)
      pass_flag <- mean_cor > 0.90

      cat(sprintf("  Mean per-taxon rho:           %.4f %s\n",
                  mean_cor, if (pass_flag) "[PASS]" else "[MODERATE]"))
      cat(sprintf("  Proportion taxa with rho≥0.90: %.1f%%\n",
                  100 * prop_pass))

      add_result(taxa_level, "4", "Abundance Concordance",
                 "Mean per-taxon Spearman rho", mean_cor, "> 0.90", pass_flag)
      add_result(taxa_level, "4", "Abundance Concordance",
                 "Proportion taxa rho >= 0.90", prop_pass, "> 0.80",
                 prop_pass > 0.80)

      # Full results table
      cor_tbl <- tibble(
        Taxa  = common_taxa,
        Spearman_rho = taxon_cors,
        Pass  = taxon_cors >= 0.90
      ) %>% arrange(Spearman_rho)

      write_tsv(cor_tbl,
                file.path(level_out, sprintf("tier4_%s_abundance_concordance.tsv",
                                             tolower(taxa_level))))

      # Histogram
      p_cor <- ggplot(cor_tbl, aes(x = Spearman_rho)) +
        geom_histogram(bins = 40, fill = "steelblue", color = "white", alpha = 0.8) +
        geom_vline(xintercept = mean_cor, color = "red",
                   linetype = "dashed", linewidth = 0.9) +
        geom_vline(xintercept = 0.90, color = "darkgreen",
                   linetype = "dotted", linewidth = 0.9) +
        annotate("text", x = mean_cor, y = Inf,
                 label = sprintf("Mean rho = %.3f", mean_cor),
                 hjust = -0.1, vjust = 1.5, color = "red", size = 3.5) +
        annotate("text", x = 0.90, y = Inf,
                 label = "rho = 0.90",
                 hjust = 1.1, vjust = 1.5, color = "darkgreen", size = 3.5) +
        labs(title = sprintf("Per-taxon abundance concordance (%s)", taxa_level),
             x = "Spearman rho (before vs after)", y = "Taxa") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      ggsave(file.path(level_plot, sprintf("17_%s_abundance_concordance.pdf",
                                           tolower(taxa_level))),
             p_cor, width = 8, height = 5)
      cat("  ✓ Abundance concordance plot saved\n")

      # Top 20 discordant taxa (lowest rho) for diagnostic use
      if (any(cor_tbl$Spearman_rho < 0.90, na.rm = TRUE)) {
        cat(sprintf("\n  Top 20 least-concordant %s:\n", tolower(taxa_level)))
        print(head(cor_tbl, 20), n = 20)
      }

    } else {
      cat(sprintf("  ⚠ Insufficient overlap (%d taxa, %d samples) — skipping\n",
                  length(common_taxa), length(common_samples)))
    }
  }, error = function(e) cat("  ⚠ Abundance concordance failed:", conditionMessage(e), "\n"))

  cat(sprintf("\n✅ All tiers complete for %s level\n", taxa_level))

}  # end taxa_level loop

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n================================================================================\n")
cat("VALIDATION SUMMARY\n")
cat("================================================================================\n\n")

write_csv(summary_results, file.path(OUTPUT_DIR, "validation_summary.csv"))

for (lv in c("ALL", TAXONOMIC_LEVELS)) {
  res_lv <- summary_results %>% filter(Level == lv)
  if (nrow(res_lv) == 0) next
  cat(sprintf("\n### %s ###\n", lv))
  print(res_lv %>% select(Tier, Analysis, Metric, Value, Criteria, Status), n = Inf)
  n_pass <- sum(res_lv$Status == "PASS", na.rm = TRUE)
  cat(sprintf("Pass rate: %d/%d (%.1f%%)\n",
              n_pass, nrow(res_lv), 100 * n_pass / nrow(res_lv)))
}

n_pass_total <- sum(summary_results$Status == "PASS", na.rm = TRUE)
n_total      <- nrow(summary_results)
cat(sprintf("\nOverall: %d/%d tests passed (%.1f%%)\n",
            n_pass_total, n_total, 100 * n_pass_total / max(n_total, 1)))

cat("\n✅ VALIDATION COMPLETE\n")
cat("Results:", OUTPUT_DIR, "\n")
cat("Figures:", PLOTS_DIR, "\n\n")
