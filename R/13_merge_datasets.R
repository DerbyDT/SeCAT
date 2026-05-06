#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   13_merge_datasets.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    13 - Dataset Merging (MetaASV Construction)
# PURPOSE:  Merge trimmed ASV sequences across studies into a unified MetaASV dataset
#
# OVERVIEW:
#   This is the final harmonisation step of the SeCAT pipeline. It loads trimmed
#   (standardised) FASTA sequences from all successful studies, deduplicates
#   identical sequences into MetaASVs, builds a cross-study feature table with
#   correct abundance accounting, assigns taxonomy using a length-weighted
#   confidence scoring system, and merges metadata across heterogeneous studies.
#   Both trimmed (post-consensus) and untrimmed (pre-consensus) datasets are
#   produced to enable downstream comparison of harmonisation effects.
#
# INPUTS:
#   - output/standardized_datasets/{study}_standardized.fasta
#       Trimmed FASTA files per study (from Stage 12)
#   - output/standardized_datasets/trim_summary.csv
#       Trimming outcome per study
#   - secat_manifest_clean.tsv
#       Manifest linking study names to feature table, taxonomy, metadata, and
#       FASTA file paths
#   - Per-study feature tables, taxonomy files, and metadata files
#       (paths specified in the manifest)
#
# OUTPUTS:
#   - output/meta_analysis/combined_sequences.fasta
#       Deduplicated MetaASV FASTA (unique sequences only)
#   - output/meta_analysis/combined_feature_table.tsv
#       Merged feature table (MetaASV x sample) with summed abundances
#   - output/meta_analysis/combined_taxonomy.tsv
#       Length-weighted best taxonomy per MetaASV
#   - output/meta_analysis/combined_metadata.tsv
#       Harmonised metadata across all studies
#   - output/meta_analysis/asv_mapping_final.tsv
#       Full mapping from original ASV IDs to MetaASV IDs
#   - output/comparison/post_consensus/  (trimmed outputs)
#   - output/comparison/pre_consensus/   (untrimmed originals for comparison)
#
# DEPENDENCIES:
#   - dplyr, tidyr, ggplot2, readr, stringr, purrr, tibble, here
#   - Biostrings (sequence I/O and manipulation)
#   - data.table (robust table parsing)
#   - R/secat_config.R (pipeline configuration)
#
# CALLED BY:
#   - Nextflow MERGE_DATASETS module
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
  library(Biostrings)
  library(data.table)
})

cat("\n=== MetaASV DATASET BUILDER ===\n")

# --- Load Pipeline Configuration ---
source(here("R/secat_config.R"))
if (!exists("SECAT_MANIFEST_PATH")) SECAT_MANIFEST_PATH <- "secat_manifest_clean.tsv"

# Print config summary
cat("=============================================================\n")
cat("  MESAP Pipeline Configuration Summary\n")
cat("=============================================================\n")
cat(sprintf("  Reference Database:      %s\n", basename(REFERENCE_DB_PATH)))
cat(sprintf("  Analysis Mode:           %s\n", ANALYSIS_MODE))
cat(sprintf("  Align ALL ASVs?          %s\n", USE_ALL_ASVS))
cat(sprintf("  Changepoint Method:      %s\n", CHANGEPOINT_PENALTY_METHOD))
cat(sprintf("  Alignment Method:        %s\n", STUDY_ALIGNMENT_METHOD))
cat(sprintf("  Simulations per primer:  %d\n", NUM_SIMULATIONS_PER_PRIMER))
cat(sprintf("  Trim increment:          %.0f bp\n", TRIM_INCREMENT))   # %d -> %.0f (numeric type)
cat("=============================================================\n")
cat("  -> Simulation mode: Custom (Grinder-inspired)\n")
cat(sprintf("     PCR bias: %s | Errors: %s | Chimeras: %s\n",
            SIMULATION_ADD_PCR_BIAS,
            SIMULATION_ADD_SEQUENCING_ERRORS,
            SIMULATION_ADD_CHIMERAS))
cat("✓ Config\n\n")

# --- Set Up Directory Paths ---
OUTPUT_DIR <- file.path(OUTDIR, "meta_analysis")
STD_DIR <- file.path(OUTDIR, "standardized_datasets")
COMPARISON_DIR <- file.path(OUTDIR, "comparison")

for (d in c(OUTPUT_DIR, COMPARISON_DIR,
            file.path(COMPARISON_DIR, "post_consensus"),
            file.path(COMPARISON_DIR, "pre_consensus"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Identify Successful Studies ---
# Determine which studies to merge. In roster mode the roster is trusted
# directly; otherwise filter by trim_summary SUCCESS status.
SELECTION_MODE <- Sys.getenv("SECAT_SELECTION_MODE", "file")
if (SELECTION_MODE == "roster") {
  roster_path <- Sys.getenv("SECAT_SELECTION_FILE", "")
  if (!file.exists(roster_path)) stop("FATAL: SECAT_SELECTION_FILE not found: ", roster_path)
  selected_studies <- readLines(roster_path)
  selected_studies <- selected_studies[!grepl("^#", selected_studies) & nchar(trimws(selected_studies)) > 0]
  cat(sprintf("Roster-selected %d studies from: %s\n", length(selected_studies), roster_path))
} else {
  selection_file <- file.path(AGGREGATED_DATA_DIR, "selected_studies_for_trim.txt")
  if (!file.exists(selection_file)) stop("FATAL: Selection file not found: ", selection_file)
  selected_studies <- readLines(selection_file)
  cat(sprintf("Loaded %d selected studies from file\n", length(selected_studies)))
}
manifest <- read_tsv(SECAT_MANIFEST_PATH, show_col_types = FALSE)

trim_summary_file <- file.path(STD_DIR, "trim_summary.csv")
trim_summary <- read_csv(trim_summary_file, show_col_types = FALSE)

if (SELECTION_MODE == "roster") {
  # Roster mode: trust the roster entirely, don't filter by trim_summary status
  successful_studies <- selected_studies
  failed_studies     <- character(0)
  cat(sprintf("Roster mode: trusting roster directly (%d studies).\n",
              length(successful_studies)))
} else {
  successful_studies <- trim_summary %>%
    filter(status == "SUCCESS") %>%
    pull(study_name)
  failed_studies <- trim_summary %>%
    filter(status != "SUCCESS") %>%
    pull(study_name)
}

cat(sprintf("Selected: %d | Passed: %d | Failed: %d\n",
            length(selected_studies),
            length(successful_studies),
            length(failed_studies)))

if (length(failed_studies) > 0) {
  cat(sprintf("  Excluded: %s\n", paste(failed_studies, collapse = ", ")))
}
cat("\n")

if (length(successful_studies) == 0) {
  stop("FATAL: No studies passed trimming!")
}

# ==============================================================================
# UTILITY FUNCTION: BioSample Run Fetching
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: fetch_run_to_biosample
# Purpose:  Map SRA run accessions to BioSample IDs via NCBI E-utilities
#
# Parameters:
#   @param run_ids   [character] - Vector of SRA run accessions (SRR/DRR/ERR)
#   @param cache_path [character] - Path to cache file for storing/loading results
#
# Returns:
#   @return [tibble] - Two-column table: run_accession, SampleID (BioSample)
#
# Notes:
#   Uses NCBI efetch in chunks of 200 to respect rate limits. Results are cached
#   to disk to avoid redundant API calls on re-runs. This is needed when feature
#   table column headers are run accessions rather than BioSample IDs, which
#   prevents direct matching against the metadata table.
# ------------------------------------------------------------------------------
fetch_run_to_biosample <- function(run_ids, cache_path) {
  if (file.exists(cache_path)) {
    cat("  ✓ Using cached run→biosample map:", cache_path, "\n")
    return(read_tsv(cache_path, show_col_types = FALSE))
  }

  cat("  Fetching run→biosample mapping from NCBI for", length(run_ids), "run IDs...\n")

  chunk_size <- 200
  chunks     <- split(run_ids, ceiling(seq_along(run_ids) / chunk_size))
  results    <- list()

  for (i in seq_along(chunks)) {
    ids_str <- paste(chunks[[i]], collapse = ",")
    url     <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
      "?db=sra&id=", ids_str, "&rettype=runinfo&retmode=text"
    )
    tryCatch({
      raw        <- read_csv(url, show_col_types = FALSE)
      results[[i]] <- raw %>%
        select(run_accession = Run, SampleID = BioSample) %>%
        filter(!is.na(run_accession), run_accession != "")
      Sys.sleep(0.4)  # respect NCBI rate limits
    }, error = function(e) {
      cat("  ✗ efetch failed for chunk", i, ":", conditionMessage(e), "\n")
    })
  }

  mapping <- bind_rows(results)
  if (nrow(mapping) > 0) {
    write_tsv(mapping, cache_path)
    cat(sprintf("  ✓ Cached %d mappings to %s\n", nrow(mapping), cache_path))
  }
  return(mapping)
}


# ==============================================================================
# UTILITY FUNCTION: ROBUST TABLE LOADING
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: load_flexible_table
# Purpose:  Robustly load TSV feature tables, handling QIIME2 biom-derived
#           headers and other format quirks
#
# Parameters:
#   @param filepath       [character] - Path to TSV file
#   @param first_col_name [character] - Name to assign to the first column
#                                       (default "ASV_ID")
#   @param as_character   [logical]   - If TRUE, read all columns as character
#
# Returns:
#   @return [tibble] - Cleaned table with standardised column names
#
# Notes:
#   Handles QIIME2 "# Constructed from biom" headers by detecting and skipping
#   them. Uses data.table::fread for speed and robustness with malformed files.
#   Strips leading "# " from column names and removes empty rows/columns.
# ------------------------------------------------------------------------------
load_flexible_table <- function(filepath, first_col_name = "ASV_ID", as_character = FALSE) {
  if (!file.exists(filepath)) return(NULL)

  # Read first few lines to detect format
  lines <- readLines(filepath, n = 5)

  # Detect QIIME2 format (biom-derived tables have a comment header)
  skip <- 0
  is_qiime2 <- FALSE

  for (i in seq_along(lines)) {
    if (grepl("^# Constructed from biom", lines[i])) {
      is_qiime2 <- TRUE
      skip <- i
      break
    }
  }

  # Use data.table::fread for robust reading
  dt <- suppressWarnings(
    fread(
      filepath,
      skip = skip,
      check.names = FALSE,
      data.table = FALSE,
      sep = "\t",
      fill = TRUE,
      quote = "",
      strip.white = TRUE,
      blank.lines.skip = TRUE,
      colClasses = if (as_character) "character" else NULL
    )
  )

  # Clean column names
  colnames(dt) <- gsub("^#\\s*", "", colnames(dt))
  colnames(dt) <- trimws(colnames(dt))

  # Set first column name
  colnames(dt)[1] <- first_col_name

  # Remove completely empty columns
  empty_cols <- sapply(dt, function(x) all(is.na(x) | x == "" | x == " "))
  if (any(empty_cols)) {
    dt <- dt[, !empty_cols, drop = FALSE]
  }

  # Remove rows with empty/NA first column
  if (nrow(dt) > 0) {
    dt <- dt[!is.na(dt[[1]]) & dt[[1]] != "" & dt[[1]] != " ", ]
  }

  return(as_tibble(dt))
}

# ------------------------------------------------------------------------------
# Function: load_taxonomy
# Purpose:  Load taxonomy file, handling multiple format variants
#
# Parameters:
#   @param filepath [character] - Path to taxonomy TSV
#
# Returns:
#   @return [tibble] - Standardised three-column table: ASV_ID, Taxon, Confidence
#
# Notes:
#   Supports three formats:
#     Case 1: Standard QIIME2 (ASV_ID + Taxon + Confidence)
#     Case 2: Legacy MetaASV format (Meta_ID used as ASV_ID)
#     Case 3: Reordered columns with same names
# ------------------------------------------------------------------------------
load_taxonomy <- function(filepath) {
  tax <- read_tsv(filepath, show_col_types = FALSE)
  cols <- colnames(tax)

  # Case 1: Standard format (ASV_ID + Taxon)
  if (all(c("ASV_ID", "Taxon") %in% cols)) {
    tax_std <- tax %>% select(ASV_ID, Taxon, Confidence)

  # Case 2: Legacy MetaASV (Taxon string, Confidence, MetaASV_ID)
  } else if (all(c("ASV_ID", "Confidence", "Meta_ID") %in% cols)) {
    tax_std <- tax %>%
      transmute(
        ASV_ID = Meta_ID,
        Taxon  = ASV_ID,
        Confidence = Confidence
      )

  # Case 3: Legacy with different column order
  } else if (all(c("Taxon", "Confidence", "ASV_ID") %in% cols)) {
    tax_std <- tax %>%
      transmute(
        ASV_ID = ASV_ID,
        Taxon  = Taxon,
        Confidence = Confidence
      )

  } else {
    stop("Unrecognized taxonomy format")
  }

  return(tax_std)
}

# ------------------------------------------------------------------------------
# Function: load_metadata
# Purpose:  Load sample metadata with all columns as character type
#
# Parameters:
#   @param filepath [character] - Path to metadata TSV
#
# Returns:
#   @return [tibble] - Metadata table with SampleID as first column
# ------------------------------------------------------------------------------
load_metadata <- function(filepath) {
  meta <- load_flexible_table(filepath, "SampleID", as_character = TRUE)
  return(meta)
}

# ==============================================================================
# DATA LOADING
# ==============================================================================
# Load feature tables, taxonomy, metadata, and trimmed sequences for each
# successful study. Sample IDs are prefixed with the study name to ensure
# uniqueness across the merged dataset.
cat("--- LOADING ---\n")
trimmed_sequences <- DNAStringSet()
untrimmed_features <- list()
untrimmed_taxonomy <- list()
all_metadata <- list()

# Track read counts per study for downstream validation
study_read_counts <- list()

for (study in successful_studies) {
  cat(sprintf("[%d/%d] %s\n", match(study, successful_studies), length(successful_studies), study))

  study_row <- manifest %>% filter(study_name == study)
  if (nrow(study_row) == 0) {
    cat("  [SKIP - not in manifest]\n")
    next
  }

  # Load feature table
  orig_counts <- load_flexible_table(study_row$asv_counts_path[1])
  if (!is.null(orig_counts)) {
    sample_cols <- colnames(orig_counts)[-1]
    n_archive <- sum(grepl("^(SRR|DRR|ERR|SAM)", sample_cols))

    # Calculate total reads
    total_reads <- sum(orig_counts[, -1], na.rm = TRUE)
    study_read_counts[[study]] <- total_reads

    cat(sprintf("  ✓ Feature table: %d ASVs, %d samples (%s reads)\n",
                nrow(orig_counts), length(sample_cols),
                format(total_reads, big.mark = ",")))

    untrimmed_features[[study]] <- orig_counts
  } else {
    cat("  ✗ Feature table failed to load\n")
    next
  }

  # --- Metadata Loading and ID Reconciliation ---
  # Feature table columns may use SRA run accessions (SRR/DRR/ERR) while
  # metadata uses BioSample IDs. If match rate is below 50%, attempt to
  # remap run accessions to BioSample IDs via NCBI E-utilities.
  meta <- load_metadata(study_row$metadata_path[1])
  if (!is.null(meta)) {
    cat(sprintf("  ✓ Metadata loaded: %d samples\n", nrow(meta)))

    # Check metadata-feature table match
    ft_samples <- colnames(orig_counts)[-1]
    meta_samples <- meta$SampleID

    n_match <- sum(ft_samples %in% meta_samples)
    match_pct <- 100 * n_match / length(ft_samples)

    cat(sprintf("  → Metadata↔Feature table: %d/%d matched (%.1f%%)\n",
                n_match, length(ft_samples), match_pct))

    if (match_pct < 50) {
      cat("  ⚠️  Low match (", round(match_pct), "%). Checking if FT IDs are run accessions...\n")

      is_run_acc <- grepl("^[DES]RR[0-9]+$", ft_samples)

      if (mean(is_run_acc) >= 0.8) {
        # Derive study directory from the manifest asv_counts_path
        study_path <- dirname(study_row$asv_counts_path[1])
        cache_path <- file.path(study_path, "run_to_sample.tsv")

        run_map <- fetch_run_to_biosample(ft_samples[is_run_acc], cache_path)

        if (nrow(run_map) > 0) {
          run_to_id <- setNames(run_map$SampleID, run_map$run_accession)

          # Rename feature table columns from run accessions to BioSample IDs
          old_names <- colnames(orig_counts)
          colnames(orig_counts)[old_names %in% names(run_to_id)] <-
            run_to_id[old_names[old_names %in% names(run_to_id)]]

          # Update ft_samples to reflect renamed columns
          ft_samples <- colnames(orig_counts)[-1]

          # Re-store the renamed FT so downstream processing uses it
          untrimmed_features[[study]] <- orig_counts

          match_pct <- 100 * sum(ft_samples %in% meta_samples) / length(ft_samples)
          cat(sprintf("  ✓ After remap: %.1f%% match (%d/%d samples)\n",
                      match_pct, sum(ft_samples %in% meta_samples), length(ft_samples)))
        } else {
          cat("  ✗ efetch returned no mappings. Metadata will be absent for this study.\n")
        }
      } else {
        cat("  ✗ FT IDs don't look like run accessions — manual investigation needed.\n")
      }
    }

    # Filter metadata to only samples present in the feature table
    meta_filtered <- meta %>% filter(SampleID %in% ft_samples)

    # Filter metadata to only samples in feature table
    meta_filtered <- meta %>% filter(SampleID %in% ft_samples)
    cat(sprintf("  → Filtered metadata: %d → %d samples (kept only FT matches)\n",
                nrow(meta), nrow(meta_filtered)))

    # Add study prefix to sample IDs and a StudyID column for provenance
    meta_filtered$SampleID <- paste0(study, "_", meta_filtered$SampleID)
    meta_filtered$StudyID <- study
    all_metadata[[study]] <- meta_filtered
  }

  # Taxonomy
  orig_tax <- load_flexible_table(study_row$taxonomy_path[1], "ASV_ID")
  if (!is.null(orig_tax) && "Taxon" %in% colnames(orig_tax)) {
    untrimmed_taxonomy[[study]] <- orig_tax
    cat(sprintf("  ✓ Taxonomy: %d ASVs\n", nrow(orig_tax)))
  }

  # --- Load Trimmed (Standardised) Sequences ---
  # These are the consensus-region-only degapped sequences from Stage 12.
  # Names are prefixed with study name to ensure cross-study uniqueness.
  std_fasta <- file.path(STD_DIR, paste0(study, "_standardized.fasta"))
  if (file.exists(std_fasta)) {
    seqs <- readDNAStringSet(std_fasta)

    # Store with study prefix in names
    names(seqs) <- paste0(study, "_", names(seqs))
    trimmed_sequences <- c(trimmed_sequences, seqs)

    cat(sprintf("  ✓ Trimmed sequences: %d\n", length(seqs)))
  }

  cat("\n")
}

# ==============================================================================
# READ COUNT CHECKPOINT 1: ORIGINAL TOTALS
# ==============================================================================
# Track total reads at each stage to detect inflation or excessive loss
# during the merge process.
cat("\n================================================================================\n")
cat("READ COUNT TRACKING - CHECKPOINT 1: ORIGINAL DATA\n")
cat("================================================================================\n")

total_original_reads <- sum(unlist(study_read_counts))
cat(sprintf("Total original reads: %s\n", format(total_original_reads, big.mark = ",")))

for (study in names(study_read_counts)) {
  pct <- 100 * study_read_counts[[study]] / total_original_reads
  cat(sprintf("  %s: %s (%.1f%%)\n",
              study,
              format(study_read_counts[[study]], big.mark = ","),
              pct))
}
cat("\n")

# ==============================================================================
# DEDUPLICATE TO MetaASVs
# ==============================================================================
# Exact sequence deduplication: identical trimmed sequences from different
# studies (or different ASVs within a study) are collapsed into a single
# MetaASV. This is the core cross-study harmonisation step.
cat("--- DEDUPLICATING ---\n")

cat(sprintf("  Sequences loaded: %d\n", length(trimmed_sequences)))

# Find unique sequences and assign MetaASV IDs
unique_seqs <- unique(trimmed_sequences)
meta_ids <- paste0("MetaASV_", seq_along(unique_seqs))
names(unique_seqs) <- meta_ids

cat(sprintf("  Unique sequences: %d\n", length(unique_seqs)))
cat(sprintf("  Compression: %.1fx\n", length(trimmed_sequences) / length(unique_seqs)))

# Diagnostic: quantify cross-study sequence sharing
seq_counts <- table(as.character(trimmed_sequences))
n_duplicated <- sum(seq_counts > 1)
cat(sprintf("  Sequences appearing >1 time: %d (%.1f%%)\n",
            n_duplicated, 100 * n_duplicated / length(seq_counts)))

if (n_duplicated > 0) {
  max_copies <- max(seq_counts)
  cat(sprintf("  Max copies of one sequence: %d\n", max_copies))

  # Show example
  most_duplicated <- names(seq_counts)[which.max(seq_counts)]
  cat(sprintf("  Most duplicated sequence: %d bp (appears %d times)\n",
              nchar(most_duplicated), max_copies))
}

# --- Build ASV-to-MetaASV Mapping Table ---
# Maps every original study-prefixed ASV ID to its deduplicated MetaASV ID
# via exact sequence matching.
final_map <- tibble(
  Original_ID = names(trimmed_sequences),
  Sequence = as.character(trimmed_sequences)
) %>%
  left_join(
    tibble(Sequence = as.character(unique_seqs), Meta_ID = meta_ids),
    by = "Sequence"
  )

# Persist the full mapping for provenance and downstream re-analysis
write_tsv(final_map, file.path(OUTPUT_DIR, "asv_mapping_final.tsv"))

# Verify mapping integrity
cat(sprintf("\n  Mapping table: %d rows\n", nrow(final_map)))
n_na_meta <- sum(is.na(final_map$Meta_ID))
if (n_na_meta > 0) {
  cat(sprintf("  ⚠️ WARNING: %d rows with NA Meta_ID!\n", n_na_meta))
}

n_dup_originals <- sum(duplicated(final_map$Original_ID))
if (n_dup_originals > 0) {
  cat(sprintf("  ⚠️ WARNING: %d duplicate Original_IDs!\n", n_dup_originals))
}

cat(sprintf("\n  %d → %d MetaASVs\n", length(trimmed_sequences), length(unique_seqs)))

# ==============================================================================
# BUILD TRIMMED FEATURE TABLE (WITH READ TRACKING)
# ==============================================================================
# For each study: prefix sample IDs, map ASV IDs to MetaASV IDs, then
# aggregate abundances by MetaASV (summing counts when multiple ASVs from
# the same study collapse to one MetaASV). Read counts are tracked at each
# step to detect inflation or loss.
cat("\n--- TRIMMED FEATURE TABLE ---\n")
feature_list <- list()
study_mapping_stats <- list()

for (study in names(untrimmed_features)) {
  cat(sprintf("  %s...\n", study))

  ft <- untrimmed_features[[study]]

  # Track original total
  original_total <- sum(ft[, -1], na.rm = TRUE)

  # Add study prefix to sample column names for cross-study uniqueness
  sample_cols <- setdiff(colnames(ft), "ASV_ID")
  new_sample_names <- paste0(study, "_", sample_cols)

  ft_prefixed <- ft
  colnames(ft_prefixed)[colnames(ft_prefixed) != "ASV_ID"] <- new_sample_names

  # Add study prefix to ASV_IDs to match final_map
  ft_prefixed$ASV_ID <- paste0(study, "_", ft_prefixed$ASV_ID)

  # Map to MetaASVs (inner join: only ASVs present in the trimmed set)
  ft_mapped <- merge(ft_prefixed, final_map[, c("Original_ID", "Meta_ID")],
                     by.x = "ASV_ID", by.y = "Original_ID", all.x = FALSE)
  ft_mapped <- ft_mapped[!is.na(ft_mapped$Meta_ID), ]

  if (nrow(ft_mapped) == 0) {
    cat("    SKIP - no mapped ASVs\n")
    next
  }

  # Check for duplicate mappings
  n_duplicates <- sum(duplicated(ft_mapped$ASV_ID))
  if (n_duplicates > 0) {
    cat(sprintf("    ⚠️ WARNING: %d duplicate ASV mappings detected!\n", n_duplicates))
    cat("    Removing duplicates...\n")
    ft_mapped <- ft_mapped[!duplicated(ft_mapped$ASV_ID), ]
  }

  # Remove ASV_ID column before aggregation
  ft_mapped$ASV_ID <- NULL
  numeric_cols <- names(ft_mapped)[sapply(ft_mapped, is.numeric)]

  # Aggregate by MetaASV: sum abundances when multiple ASVs map to one MetaASV
  ft_agg <- aggregate(ft_mapped[, numeric_cols, drop = FALSE],
                      by = list(Meta_ID = ft_mapped$Meta_ID),
                      FUN = sum, na.rm = TRUE)

  # Verify totals (retention should be ~100% if all ASVs were trimmed)
  mapped_total <- sum(ft_agg[, -1], na.rm = TRUE)
  retention_pct <- 100 * mapped_total / original_total

  cat(sprintf("    Original reads: %s\n", format(original_total, big.mark = ",")))
  cat(sprintf("    Mapped reads:   %s\n", format(mapped_total, big.mark = ",")))
  cat(sprintf("    Retention:      %.2f%%\n", retention_pct))
  cat(sprintf("    MetaASVs:       %d\n", nrow(ft_agg)))

  if (retention_pct > 101) {
    cat(sprintf("    ⚠️ READ INFLATION: %.2f%% > 100%%\n", retention_pct))
  } else if (retention_pct < 80) {
    cat(sprintf("    ⚠️ LOW RETENTION: %.2f%% < 80%%\n", retention_pct))
  }

  # Store stats
  study_mapping_stats[[study]] <- list(
    original = original_total,
    mapped = mapped_total,
    retention = retention_pct,
    n_metaasvs = nrow(ft_agg)
  )

  feature_list[[study]] <- ft_agg
}

# ==============================================================================
# READ COUNT CHECKPOINT 2: AFTER MAPPING
# ==============================================================================
cat("\n================================================================================\n")
cat("READ COUNT TRACKING - CHECKPOINT 2: AFTER MetaASV MAPPING\n")
cat("================================================================================\n")

total_mapped <- sum(sapply(study_mapping_stats, function(x) x$mapped))
overall_retention <- 100 * total_mapped / total_original_reads

cat(sprintf("Total original reads: %s\n", format(total_original_reads, big.mark = ",")))
cat(sprintf("Total mapped reads:   %s\n", format(total_mapped, big.mark = ",")))
cat(sprintf("Overall retention:    %.2f%%\n\n", overall_retention))

cat("Per-study breakdown:\n")
for (study in names(study_mapping_stats)) {
  stats <- study_mapping_stats[[study]]
  cat(sprintf("  %s:\n", study))
  cat(sprintf("    Original: %s\n", format(stats$original, big.mark = ",")))
  cat(sprintf("    Mapped:   %s (%.2f%%)\n", format(stats$mapped, big.mark = ","), stats$retention))
  cat(sprintf("    MetaASVs: %d\n", stats$n_metaasvs))
}
cat("\n")

if (overall_retention > 100.5) {
  cat("⚠️ WARNING: Read count inflation detected (>100.5%)!\n")
  cat("   This suggests duplicate counting during aggregation.\n\n")
} else if (overall_retention < 80) {
  cat("⚠️ WARNING: Low read retention (<80%)!\n")
  cat("   Many ASVs may have been filtered during trimming.\n\n")
}

# ==============================================================================
# MERGE FEATURE TABLES
# ==============================================================================
# Sequentially full-outer-join per-study feature tables by Meta_ID so that
# all MetaASVs and all samples are represented. Missing values (a MetaASV
# not observed in a given study) are filled with zero.
cat("--- MERGING FEATURE TABLES ---\n")
cat(sprintf("  Merging %d studies: ", length(feature_list)))

feature_merged <- feature_list[[1]]

if (length(feature_list) > 1) {
  for (i in 2:length(feature_list)) {
    cat(sprintf("%d ", i))
    feature_merged <- merge(feature_merged, feature_list[[i]], by = "Meta_ID", all = TRUE)
  }
}
cat("✓\n")

# Replace NAs with 0 (absent MetaASVs in a study have zero counts)
for (col in names(feature_merged)) {
  if (is.numeric(feature_merged[[col]])) {
    feature_merged[[col]][is.na(feature_merged[[col]])] <- 0
  }
}

trimmed_feature_final <- feature_merged
names(trimmed_feature_final)[names(trimmed_feature_final) == "Meta_ID"] <- "ASV_ID"

# ==============================================================================
# READ COUNT CHECKPOINT 3: FINAL MERGED TABLE
# ==============================================================================
cat("\n================================================================================\n")
cat("READ COUNT TRACKING - CHECKPOINT 3: FINAL MERGED TABLE\n")
cat("================================================================================\n")

final_total <- sum(trimmed_feature_final[, -1], na.rm = TRUE)
final_retention <- 100 * final_total / total_original_reads

cat(sprintf("Total original reads:  %s\n", format(total_original_reads, big.mark = ",")))
cat(sprintf("Total final reads:     %s\n", format(final_total, big.mark = ",")))
cat(sprintf("Final retention:       %.2f%%\n", final_retention))
cat(sprintf("Final MetaASVs:        %d\n", nrow(trimmed_feature_final)))
cat(sprintf("Final samples:         %d\n", ncol(trimmed_feature_final) - 1))
cat("\n")

if (final_retention > 100.1) {
  cat("❌ VALIDATION FAILED: Read inflation detected!\n")
  cat("   Expected: ≤100%, Observed: %.2f%%\n", final_retention)
  cat("   Difference: +%s reads\n", format(final_total - total_original_reads, big.mark = ","))
} else if (final_retention >= 95 && final_retention <= 100) {
  cat("✅ VALIDATION PASSED: Excellent read retention (%.2f%%)\n", final_retention)
} else if (final_retention >= 80 && final_retention < 95) {
  cat("⚠️ WARNING: Moderate read loss (%.2f%%)\n", final_retention)
  cat("   Lost: %s reads\n", format(total_original_reads - final_total, big.mark = ","))
} else {
  cat("❌ VALIDATION FAILED: Excessive read loss (%.2f%%)\n", final_retention)
  cat("   Lost: %s reads\n", format(total_original_reads - final_total, big.mark = ","))
}
cat("\n")

# ==============================================================================
# TAXONOMY WITH LENGTH-WEIGHTED SELECTION
# ==============================================================================
# When multiple ASVs from different studies collapse to the same MetaASV,
# we need to pick the "best" taxonomy. This section implements a composite
# scoring system: 70% classifier confidence + 30% original (pre-trim)
# sequence length. Longer source sequences tend to yield more accurate
# taxonomic classifications because the classifier has more informative
# positions to work with.
cat("--- TAXONOMY (Length-Weighted) ---\n")

# --- Step 1: Load Original Sequence Lengths ---
# We use the ORIGINAL (pre-trim) sequence lengths, not trimmed lengths,
# because the original length reflects how much taxonomic information
# the classifier had when assigning the taxonomy.
cat("  Step 1: Loading original sequence lengths...\n")

original_lengths <- list()

for (study in names(untrimmed_features)) {
  study_row <- manifest %>% filter(study_name == study)
  if (nrow(study_row) == 0) next

  # Load original (untrimmed) sequences
  original_fasta_path <- study_row$asv_fasta_path[1]

  if (!is.na(original_fasta_path) && file.exists(original_fasta_path)) {
    orig_seqs <- readDNAStringSet(original_fasta_path)

    # Store lengths with study prefix
    lengths_df <- tibble(
      ASV_ID = paste0(study, "_", names(orig_seqs)),
      Original_Length = width(orig_seqs)
    )

    original_lengths[[study]] <- lengths_df

    cat(sprintf("    %s: %d sequences (%.0f bp mean)\n",
                study, nrow(lengths_df), mean(lengths_df$Original_Length)))
  } else {
    cat(sprintf("    %s: ⚠️ Original sequences not found at %s\n",
                study, original_fasta_path))
  }
}

# Check if we got any lengths
if (length(original_lengths) == 0) {
  cat("\n  ⚠️ WARNING: No original sequence files found!\n")
  cat("     Falling back to using trimmed sequence lengths...\n\n")

  # Fallback: use trimmed sequences (less ideal but functional)
  all_lengths <- tibble(
    ASV_ID = names(trimmed_sequences),
    Original_Length = width(trimmed_sequences)
  )

  cat(sprintf("  Using trimmed lengths: %d sequences\n", nrow(all_lengths)))

} else {
  # Combine all length data
  all_lengths <- bind_rows(original_lengths)

  cat(sprintf("\n  Combined: %d ASVs with original lengths\n", nrow(all_lengths)))
  cat(sprintf("  Length range: %d - %d bp (mean: %.0f bp)\n",
              min(all_lengths$Original_Length),
              max(all_lengths$Original_Length),
              mean(all_lengths$Original_Length)))
}

# --- Step 2: Prepare Taxonomy Tables ---
# Add study prefix to ASV IDs so they match the mapping table
cat("\n  Step 2: Preparing taxonomy tables...\n")
all_tax_list <- list()

for (study in names(untrimmed_taxonomy)) {
  tax <- untrimmed_taxonomy[[study]]

  # Add study prefix
  tax$ASV_ID <- paste0(study, "_", tax$ASV_ID)
  tax$Study <- study

  all_tax_list[[study]] <- tax

  cat(sprintf("    %s: %d taxonomy entries\n", study, nrow(tax)))
}

# --- Step 3: Combine All Taxonomies ---
all_tax <- bind_rows(all_tax_list) %>%
  mutate(Confidence = as.numeric(Confidence))

cat(sprintf("\n  Combined taxonomy: %d entries\n", nrow(all_tax)))

# --- Step 4: Add Original Lengths to Taxonomy ---
cat("\n  Step 3: Merging lengths with taxonomy...\n")

all_tax_with_length <- all_tax %>%
  left_join(all_lengths, by = "ASV_ID")

n_with_length <- sum(!is.na(all_tax_with_length$Original_Length))
cat(sprintf("    %d/%d taxonomy entries have length data (%.1f%%)\n",
            n_with_length, nrow(all_tax_with_length),
            100 * n_with_length / nrow(all_tax_with_length)))

# For entries without length data, use median as fallback
if (n_with_length > 0) {
  median_length <- median(all_tax_with_length$Original_Length, na.rm = TRUE)
  n_missing <- sum(is.na(all_tax_with_length$Original_Length))

  if (n_missing > 0) {
    all_tax_with_length$Original_Length[is.na(all_tax_with_length$Original_Length)] <- median_length
    cat(sprintf("    Using median length (%d bp) for %d missing values\n",
                median_length, n_missing))
  }
} else {
  stop("CRITICAL: No sequence lengths available!")
}

# --- Step 5: Map Taxonomy to MetaASVs ---
cat("\n  Step 4: Mapping to MetaASVs...\n")

trimmed_tax_with_length <- all_tax_with_length %>%
  left_join(
    final_map %>% select(Original_ID, Meta_ID),
    by = c("ASV_ID" = "Original_ID"),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(Meta_ID))

cat(sprintf("    After join: %d rows\n", nrow(trimmed_tax_with_length)))

# --- Step 6: Calculate Composite Score ---
# Composite = 0.70 * normalised_confidence + 0.30 * normalised_length
# Length is normalised within each MetaASV's contributing ASVs (0-1 scale).
# The 70/30 weighting prioritises classifier confidence while giving a
# meaningful bonus to ASVs from longer amplicons.
cat("\n  Step 5: Calculating length-weighted confidence scores...\n")

trimmed_tax_scored <- trimmed_tax_with_length %>%
  group_by(Meta_ID) %>%
  mutate(
    # Normalize length to 0-1 scale (within this MetaASV's contributors)
    Length_Score = (Original_Length - min(Original_Length)) /
                   (max(Original_Length) - min(Original_Length) + 0.001),

    # Normalize confidence to 0-1 scale
    Confidence_Score = Confidence / 100,

    # Composite score: weighted average (70% confidence, 30% length)
    Composite_Score = (0.70 * Confidence_Score) + (0.30 * Length_Score)
  ) %>%
  ungroup()

# Show scoring example (only if there are MetaASVs with multiple sources)
metaasvs_with_multiple <- trimmed_tax_scored %>%
  group_by(Meta_ID) %>%
  filter(n() > 1) %>%
  pull(Meta_ID) %>%
  unique()

if (length(metaasvs_with_multiple) > 0) {
  example_metaasv <- metaasvs_with_multiple[1]

  cat("\n  Example scoring for MetaASV with multiple source sequences:\n")
  example_data <- trimmed_tax_scored %>%
    filter(Meta_ID == example_metaasv) %>%
    select(Study, Original_Length, Confidence, Length_Score,
           Confidence_Score, Composite_Score) %>%
    arrange(desc(Composite_Score))

  print(as.data.frame(example_data), row.names = FALSE)
} else {
  cat("\n  ℹ️ All MetaASVs have single source (perfect 1:1 mapping)\n")
  cat("     Length weighting will not affect results.\n")
}

# --- Step 7: Select Best Taxonomy per MetaASV ---
# For each MetaASV, pick the contributing ASV with the highest composite
# score. Ties are broken by confidence, then by original sequence length.
cat("\n  Step 6: Selecting best taxonomy per MetaASV...\n")

trimmed_tax_final <- trimmed_tax_scored %>%
  group_by(Meta_ID) %>%
  arrange(desc(Composite_Score), desc(Confidence), desc(Original_Length)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    ASV_ID = Meta_ID,
    Taxon,
    Confidence,
    Original_Length,
    Composite_Score,
    Source_Study = Study
  )

cat(sprintf("\n  ✓ Final: %d/%d MetaASVs with taxonomy (%.1f%%)\n",
            nrow(trimmed_tax_final),
            length(unique_seqs),
            100 * nrow(trimmed_tax_final) / length(unique_seqs)))

# --- Step 8: Length-Weighting Diagnostics ---
# Compare composite-score selection against confidence-only selection to
# quantify how much the length weighting changed taxonomy assignments.
if (length(metaasvs_with_multiple) > 0) {
  cat("\n  Step 7: Length-weighting diagnostics...\n")

  # Compare what would have been selected by confidence alone
  confidence_only <- trimmed_tax_scored %>%
    group_by(Meta_ID) %>%
    arrange(desc(Confidence)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(Meta_ID, Conf_Study = Study, Conf_Length = Original_Length)

  composite_selection <- trimmed_tax_final %>%
    select(Meta_ID = ASV_ID, Comp_Study = Source_Study, Comp_Length = Original_Length)

  comparison <- confidence_only %>%
    inner_join(composite_selection, by = "Meta_ID")

  n_different <- sum(comparison$Conf_Study != comparison$Comp_Study)
  n_longer <- sum(comparison$Comp_Length > comparison$Conf_Length, na.rm = TRUE)

  cat(sprintf("    MetaASVs where selection changed: %d/%d (%.1f%%)\n",
              n_different, nrow(comparison),
              100 * n_different / nrow(comparison)))

  if (n_different > 0) {
    cat(sprintf("    MetaASVs now using longer sequence: %d/%d (%.1f%%)\n",
                n_longer, nrow(comparison),
                100 * n_longer / nrow(comparison)))

    cat("\n    Length distribution of changes:\n")

    changed <- comparison %>%
      filter(Conf_Study != Comp_Study) %>%
      mutate(Length_Gain = Comp_Length - Conf_Length)

    cat(sprintf("      Mean length gain: %.0f bp\n", mean(changed$Length_Gain, na.rm = TRUE)))
    cat(sprintf("      Median length gain: %.0f bp\n", median(changed$Length_Gain, na.rm = TRUE)))
    cat(sprintf("      Max length gain: %.0f bp\n", max(changed$Length_Gain, na.rm = TRUE)))

    # Show example of improved selection
    example_improvement <- changed %>%
      arrange(desc(Length_Gain)) %>%
      slice_head(n = 1)

    if (nrow(example_improvement) > 0) {
      cat("\n    Example improvement:\n")
      cat(sprintf("      Old: %s (%d bp)\n",
                  example_improvement$Conf_Study,
                  example_improvement$Conf_Length))
      cat(sprintf("      New: %s (%d bp) [+%d bp]\n",
                  example_improvement$Comp_Study,
                  example_improvement$Comp_Length,
                  example_improvement$Length_Gain))
    }
  }
} else {
  cat("\n  Step 7: Skipping diagnostics (no competing taxonomies)\n")
}

# --- Step 9: Validation Checks ---
if (nrow(trimmed_tax_final) == 0) {
  cat("  ❌ CRITICAL: No taxonomy mapped!\n")
  stop("Taxonomy mapping completely failed!")
}

if (nrow(trimmed_tax_final) < length(unique_seqs) * 0.90) {
  cat(sprintf("  ⚠️ WARNING: Only %.1f%% coverage (expected >90%%)\n",
              100 * nrow(trimmed_tax_final) / length(unique_seqs)))
}

cat("\n✅ Length-weighted taxonomy assignment complete!\n\n")

# ==============================================================================
# UNTRIMMED DATASET
# ==============================================================================
# Build a parallel untrimmed (pre-consensus) feature table for comparison.
# This allows downstream analyses to quantify the effect of harmonisation
# on community composition metrics.
cat("\n--- UNTRIMMED DATASET ---\n")

untrimmed_prefixed <- list()
untrimmed_read_totals <- list()

for (study in names(untrimmed_features)) {
  ft <- untrimmed_features[[study]]
  sample_cols <- setdiff(colnames(ft), "ASV_ID")
  new_sample_names <- paste0(study, "_", sample_cols)

  ft_renamed <- ft
  colnames(ft_renamed)[colnames(ft_renamed) != "ASV_ID"] <- new_sample_names
  untrimmed_prefixed[[study]] <- ft_renamed

  # Track reads
  total_reads <- sum(ft_renamed[, -1], na.rm = TRUE)
  untrimmed_read_totals[[study]] <- total_reads

  cat(sprintf("  %s: %d samples, %s reads\n",
              study, length(sample_cols), format(total_reads, big.mark = ",")))
}

# Full outer join of untrimmed feature tables across studies
untrimmed_feature_final <- untrimmed_prefixed[[1]]
if (length(untrimmed_prefixed) > 1) {
  for (i in 2:length(untrimmed_prefixed)) {
    untrimmed_feature_final <- merge(untrimmed_feature_final,
                                      untrimmed_prefixed[[i]],
                                      by = "ASV_ID", all = TRUE)
  }
}

for (col in names(untrimmed_feature_final)) {
  if (is.numeric(untrimmed_feature_final[[col]])) {
    untrimmed_feature_final[[col]][is.na(untrimmed_feature_final[[col]])] <- 0
  }
}

# Combine taxonomy from all studies (deduplicate by ASV_ID)
untrimmed_tax_final <- bind_rows(untrimmed_taxonomy) %>%
  distinct(ASV_ID, .keep_all = TRUE)

untrimmed_total_reads <- sum(untrimmed_feature_final[, -1], na.rm = TRUE)

cat(sprintf("\n  ✓ %d original ASVs × %d samples\n",
            nrow(untrimmed_feature_final), ncol(untrimmed_feature_final) - 1))
cat(sprintf("  ✓ Total reads: %s\n", format(untrimmed_total_reads, big.mark = ",")))

# ==============================================================================
# SAVE ALL OUTPUTS
# ==============================================================================
# Write combined outputs for downstream analysis:
#   - meta_analysis/: the harmonised (trimmed) MetaASV dataset
#   - comparison/post_consensus/: copy of trimmed outputs for paired comparison
#   - comparison/pre_consensus/: untrimmed originals for paired comparison
cat("\n--- SAVING ---\n")

writeXStringSet(unique_seqs, file.path(OUTPUT_DIR, "combined_sequences.fasta"))
write_tsv(trimmed_feature_final, file.path(OUTPUT_DIR, "combined_feature_table.tsv"))
write_tsv(trimmed_tax_final, file.path(OUTPUT_DIR, "combined_taxonomy.tsv"))

file.copy(file.path(OUTPUT_DIR, c("combined_feature_table.tsv", "combined_taxonomy.tsv")),
          file.path(COMPARISON_DIR, "post_consensus/"), overwrite = TRUE)

write_tsv(untrimmed_feature_final,
          file.path(COMPARISON_DIR, "pre_consensus/feature_table.tsv"))
write_tsv(untrimmed_tax_final,
          file.path(COMPARISON_DIR, "pre_consensus/taxonomy.tsv"))

if (length(all_metadata) > 0) {
  # Combine metadata from all studies, deduplicating by SampleID
  master_metadata <- bind_rows(all_metadata) %>% distinct(SampleID, .keep_all = TRUE)
  write_tsv(master_metadata, file.path(OUTPUT_DIR, "combined_metadata.tsv"))
  file.copy(file.path(OUTPUT_DIR, "combined_metadata.tsv"),
            c(file.path(COMPARISON_DIR, "post_consensus/metadata.tsv"),
              file.path(COMPARISON_DIR, "pre_consensus/metadata.tsv")),
            overwrite = TRUE)

  cat(sprintf("  ✓ Metadata: %d samples\n", nrow(master_metadata)))
}

cat("  ✓ All datasets saved\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n================================================================================\n")
cat("✅ MERGE COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("  Studies included:      %d (%s)\n",
            length(successful_studies),
            paste(successful_studies, collapse = ", ")))
cat(sprintf("  MetaASVs (trimmed):    %d ASVs × %d samples\n",
            nrow(trimmed_feature_final), ncol(trimmed_feature_final) - 1))
cat(sprintf("  Original (untrimmed):  %d ASVs × %d samples\n",
            nrow(untrimmed_feature_final), ncol(untrimmed_feature_final) - 1))
cat(sprintf("  Deduplication ratio:   %.1fx\n",
            nrow(untrimmed_feature_final) / nrow(trimmed_feature_final)))

if (length(all_metadata) > 0) {
  cat(sprintf("  Metadata samples:      %d\n", nrow(master_metadata)))
}

cat("\n  Read Retention Summary:\n")
cat(sprintf("    Original:            %s reads\n", format(total_original_reads, big.mark = ",")))
cat(sprintf("    Trimmed (merged):    %s reads (%.2f%%)\n",
            format(final_total, big.mark = ","), final_retention))
cat(sprintf("    Untrimmed (merged):  %s reads\n", format(untrimmed_total_reads, big.mark = ",")))

cat("\nOutputs:\n")
cat(sprintf("  - Post-consensus (trimmed):   %s\n", file.path(COMPARISON_DIR, "post_consensus/")))
cat(sprintf("  - Pre-consensus (original): %s\n", file.path(COMPARISON_DIR, "pre_consensus/")))
cat(sprintf("  - Combined outputs:    %s\n", OUTPUT_DIR))
cat("================================================================================\n")
