
# ==============================================================================
# SCRIPT:   00_clean_data.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 0 - Data Cleaning & Quality Control
# PURPOSE:  Remove contaminant ASVs, filter metadata, and synchronise all study
#           data files prior to cross-study harmonisation.
#
# OVERVIEW:
#   This script iterates over every study listed in the SeCAT manifest and
#   performs four cleaning operations: (1) optional metadata filtering to retain
#   only samples from target environments, (2) removal of Chloroplast and
#   Mitochondria ASVs from the taxonomy table, (3) synchronisation of the
#   feature table to retain only valid ASVs and samples, and (4) subsetting
#   of the representative FASTA to match the cleaned feature table. The cleaned
#   files are written to per-study "clean/" subdirectories and an updated
#   manifest is produced for downstream pipeline stages.
#
# INPUTS:
#   - secat_manifest.tsv (or path set via SECAT_MANIFEST env var): master
#     manifest with per-study paths to ASV counts, taxonomy, metadata, and FASTA
#
# OUTPUTS:
#   - Per-study cleaned files in <study_dir>/clean/:
#       - Feature table (TSV, preserving QIIME2 headers if present)
#       - Taxonomy table (TSV)
#       - Metadata table (TSV)
#       - Representative sequences (FASTA)
#   - cleaned_data/secat_manifest_clean.tsv: updated manifest pointing to
#     cleaned files, with failed/empty studies removed
#
# DEPENDENCIES:
#   - R packages: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr,
#     forcats, Biostrings, data.table
#
# CALLED BY:
#   - Nextflow module: CLEAN_DATA (Stage 0), or run standalone
#
# v2.2 changes vs v2.1:
#   - load_feature_table() separated from load_table_robust() -- feature tables
#     require guaranteed numeric count columns; metadata/taxonomy do not.
#     Previously, colClasses="character" + a heuristic coercion loop caused
#     float-encoded counts ("3653.0") to silently remain as character in some
#     studies (notably Kardish_2023), producing all-zero colSums and dropping
#     all samples.
#   - Full per-sample read count diagnostics printed before and after filtering.
#   - QIIME2 skip logic made explicit and validated.
# ==============================================================================

# --- Load Required Packages ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(Biostrings)
  library(data.table)
})

# --- Configure I/O Paths ---
# Environment variables allow Nextflow to inject paths; fallbacks enable local runs
Read paths from environment (Nextflow mode) with fallback for direct invocation
MANIFEST_IN  <- Sys.getenv("SECAT_MANIFEST",  unset = "secat_manifest.tsv")
OUTDIR       <- Sys.getenv("SECAT_OUTDIR",     unset = ".")
MANIFEST_OUT <- file.path(OUTDIR, "cleaned_data", "secat_manifest_clean.tsv")
dir.create(file.path(OUTDIR, "cleaned_data"), recursive = TRUE, showWarnings = FALSE)

# --- Filtering Thresholds ---
# Regex for non-target organellar lineages commonly retained in 16S datasets
TAXA_TO_REMOVE       <- "Chloroplast|Mitochondria"
# Studies with fewer samples than this after filtering are excluded entirely
MIN_SAMPLES          <- 3
# Samples with fewer total reads than this are dropped (removes empties)
MIN_READS_PER_SAMPLE <- 1

message("========================================================")
message("          SeCAT DATA CLEANING UTILITY v2.2")
message("        (Robust numeric loading + diagnostics)")
message("========================================================")
message(paste("Input manifest :", MANIFEST_IN))
message(paste("Output manifest:", MANIFEST_OUT))
message(paste("Drop taxa matching:", TAXA_TO_REMOVE))
message(paste("Min samples per study:", MIN_SAMPLES))
message(paste("Min reads per sample:", MIN_READS_PER_SAMPLE))
message("========================================================")

if (!file.exists(MANIFEST_IN)) stop("FATAL: Manifest not found: ", MANIFEST_IN)

manifest       <- read_tsv(MANIFEST_IN, show_col_types = FALSE)
manifest_clean <- manifest

# ==============================================================================
# ------------------------------------------------------------------------------
# Function: load_table_robust
# Purpose:  Safely load TSV files (metadata, taxonomy) that may contain mixed
#           column types or QIIME2 comment headers.
#
# Parameters:
#   @param filepath      [character] - Path to the TSV file
#   @param first_col_name [character|NULL] - If provided, rename column 1 to this
#
# Returns:
#   @return [tibble] - Cleaned table with all columns as character type.
#     Attribute "qiime_headers" attached if QIIME2 comment lines were detected.
#
# Notes:
#   All columns are loaded as character to avoid type-inference issues with
#   mixed-type metadata columns (e.g., dates, GPS coordinates). The caller is
#   responsible for any numeric coercion. QIIME2 biom-exported files often have
#   a "# Constructed from biom file" comment line (no tabs) above the actual
#   header row -- these are skipped but preserved as an attribute for re-writing.
# ------------------------------------------------------------------------------
# ==============================================================================

load_table_robust <- function(filepath, first_col_name = NULL) {

  lines <- readLines(filepath, n = 5)

  # --- QIIME2 Comment Detection ---
  # Skip only pure QIIME2 comment lines (e.g. "# Constructed from biom file").
  # The actual column header in QIIME2 files starts with "#OTU ID" or
  # "#Feature ID" -- it contains tabs and must NOT be skipped; fread reads it
  # as the column header row. Pure comment lines have no tabs.
  skip <- 0
  qiime_headers <- NULL
  for (i in seq_along(lines)) {
    line <- lines[i]
    # Pure comments start with # but have no tab-delimited fields
    is_pure_comment <- startsWith(line, "#") && !grepl("\t", line)
    if (is_pure_comment) {
      skip <- i
      qiime_headers <- c(qiime_headers, line)
    } else {
      break
    }
  }

  # Read with all-character columns to prevent type-inference errors on mixed data
  dt <- fread(
    filepath,
    skip        = skip,
    check.names = FALSE,
    data.table  = FALSE,
    sep         = "\t",
    fill        = TRUE,
    quote       = "",
    strip.white = TRUE,
    blank.lines.skip = TRUE,
    colClasses  = "character"
  )

  # Strip leading "# " from QIIME2 header names (e.g., "#OTU ID" -> "OTU ID")
  colnames(dt) <- gsub("^#\\s*", "", colnames(dt))
  colnames(dt) <- trimws(colnames(dt))

  # Optionally standardise the first column name for downstream joins
  if (!is.null(first_col_name) && ncol(dt) >= 1) {
    colnames(dt)[1] <- first_col_name
  }

  # Remove fully empty columns (character-safe)
  empty_cols <- sapply(dt, function(x) all(is.na(x) | trimws(x) == ""))
  if (any(empty_cols)) dt <- dt[, !empty_cols, drop = FALSE]

  # Remove rows with empty/NA first column
  if (nrow(dt) > 0) {
    first <- trimws(dt[[1]])
    dt <- dt[!is.na(first) & first != "", ]
  }

  result <- as_tibble(dt)
  if (!is.null(qiime_headers)) attr(result, "qiime_headers") <- qiime_headers
  return(result)
}

# ==============================================================================
# ------------------------------------------------------------------------------
# Function: load_feature_table
# Purpose:  Load an ASV feature table with guaranteed numeric count columns.
#
# Parameters:
#   @param filepath [character] - Path to the feature table TSV
#
# Returns:
#   @return [tibble] - Feature table with ASV_ID (character) in column 1 and
#     all subsequent columns coerced to numeric. Attribute "qiime_headers"
#     attached if QIIME2 comment lines were detected.
#
# Notes:
#   Separated from load_table_robust() because feature tables have a known
#   structure: col 1 = ASV ID (character), cols 2..N = numeric counts.
#   QIIME2 biom exports sometimes encode integer counts as floats ("3653.0");
#   direct type inference would treat these correctly, but mixed-format files
#   can confuse fread. Strategy: load all as character, then explicitly coerce
#   count columns with validation. Columns where >10% of values produce NAs
#   on coercion are flagged and left as character (likely mislabelled columns).
# ------------------------------------------------------------------------------
# ==============================================================================

load_feature_table <- function(filepath) {

  lines <- readLines(filepath, n = 10)

  # --- QIIME2 Comment Detection ---
  # Skip only pure QIIME2 comment lines (no tabs -- not the header row).
  # "#OTU ID\tSRR..." is the actual header and must not be skipped.
  skip <- 0
  qiime_headers <- NULL
  for (i in seq_along(lines)) {
    line <- lines[i]
    is_pure_comment <- startsWith(line, "#") && !grepl("\t", line)
    if (is_pure_comment) {
      skip <- i
      qiime_headers <- c(qiime_headers, line)
    } else {
      break
    }
  }

  message(sprintf("       [FT loader] File: %s", basename(filepath)))
  message(sprintf("       [FT loader] Skipping %d header/comment line(s)", skip))
  if (skip > 0) message(sprintf("       [FT loader] Comment: %s", lines[1]))

  # Read with all-character to avoid any type inference crashes
  dt <- fread(
    filepath,
    skip        = skip,
    check.names = FALSE,
    data.table  = FALSE,
    sep         = "\t",
    fill        = TRUE,
    quote       = "",
    strip.white = TRUE,
    blank.lines.skip = TRUE,
    colClasses  = "character"
  )

  # Strip QIIME2 "#" prefix from column names
  colnames(dt) <- gsub("^#\\s*", "", colnames(dt))
  colnames(dt) <- trimws(colnames(dt))

  # Validate structure: need at least 2 columns (ASV ID + >=1 sample)
  if (ncol(dt) < 2) stop("Feature table has fewer than 2 columns after loading: ", filepath)

  # Standardise the ASV identifier column name
  colnames(dt)[1] <- "ASV_ID"

  # Remove fully empty rows/cols
  empty_cols <- sapply(dt, function(x) all(is.na(x) | trimws(x) == ""))
  if (any(empty_cols)) dt <- dt[, !empty_cols, drop = FALSE]
  if (nrow(dt) > 0) {
    first <- trimws(dt[[1]])
    dt <- dt[!is.na(first) & first != "", ]
  }

  # --- Numeric Coercion of Count Columns ---
  # All columns except ASV_ID must be numeric.
  # Counts may be stored as floats ("3653.0") from biom export.
  sample_cols <- names(dt)[-1]
  n_coerced   <- 0
  n_failed    <- 0

  for (col in sample_cols) {
    vals      <- dt[[col]]
    converted <- suppressWarnings(as.numeric(vals))
    na_before <- sum(is.na(vals))
    na_after  <- sum(is.na(converted))
    new_nas   <- na_after - na_before

    # Safety check: if coercion introduces >10% NAs, the column likely isn't
    # numeric data (e.g., a taxonomy column appended to the wrong file)
    if (new_nas > length(vals) * 0.1) {
      # More than 10% new NAs introduced -- coercion is mangling real data
      warning(sprintf(
        "[FT loader] Column '%s': %d new NAs after as.numeric() — keeping as character",
        col, new_nas
      ))
      n_failed <- n_failed + 1
    } else {
      dt[[col]] <- converted
      n_coerced <- n_coerced + 1
    }
  }

  message(sprintf("       [FT loader] Coerced %d/%d sample columns to numeric (%d failed)",
                  n_coerced, length(sample_cols), n_failed))

  # --- Diagnostic Summary ---
  # Report count value distribution to help catch silent loading errors
  all_counts <- unlist(dt[, -1], use.names = FALSE)
  all_counts <- suppressWarnings(as.numeric(all_counts))
  all_counts <- all_counts[!is.na(all_counts)]
  message(sprintf("       [FT loader] Count range: min=%.1f, max=%.1f, mean=%.1f",
                  min(all_counts), max(all_counts), mean(all_counts)))
  message(sprintf("       [FT loader] Zero entries: %d / %d (%.1f%%)",
                  sum(all_counts == 0), length(all_counts),
                  100 * sum(all_counts == 0) / length(all_counts)))

  result <- as_tibble(dt)
  if (!is.null(qiime_headers)) attr(result, "qiime_headers") <- qiime_headers
  return(result)
}

# ==============================================================================
# ------------------------------------------------------------------------------
# Function: make_clean_path
# Purpose:  Derive the output path for a cleaned file by inserting a "clean/"
#           subdirectory into the original file's parent directory.
#
# Parameters:
#   @param path [character] - Original file path
#
# Returns:
#   @return [character] - Path with "clean/" inserted before the filename
# ------------------------------------------------------------------------------
# ==============================================================================

make_clean_path <- function(path) {
  path       <- path[1]
  parent_dir <- dirname(path)
  fname      <- basename(path)
  file.path(parent_dir, "clean", fname)
}

# ==============================================================================
# ------------------------------------------------------------------------------
# Function: auto_detect_and_map_samples
# Purpose:  Reconcile sample ID mismatches between feature tables and metadata,
#           commonly caused by BioSample (SAMD*) vs Run accession (SRR/DRR/ERR)
#           identifier discrepancies in public repository datasets.
#
# Parameters:
#   @param ft_sample_cols [character] - Sample IDs from the feature table columns
#   @param meta_clean     [tibble]   - Filtered metadata table
#   @param sample_id_col  [character] - Name of the sample ID column in metadata
#   @param orig_metadata  [character] - Path to the original (unfiltered) metadata
#   @param study_name     [character] - Study identifier for logging
#
# Returns:
#   @return [list] with elements:
#     - valid_samples   [character] - Reconciled sample IDs matching FT columns
#     - meta_clean      [tibble]    - Metadata with sample IDs remapped if needed
#     - mapping_applied [logical]   - Whether a mapping transformation was applied
#
# Notes:
#   Uses a three-tier strategy: (1) column name heuristic to find a run
#   accession column, (2) column content scan for SRA-pattern values, (3)
#   sequential positional mapping as a last resort (with warning). The 90%
#   match threshold avoids unnecessary remapping when IDs already align.
# ------------------------------------------------------------------------------
# ==============================================================================

auto_detect_and_map_samples <- function(ft_sample_cols, meta_clean, sample_id_col,
                                        orig_metadata, study_name) {

  valid_samples <- meta_clean[[sample_id_col]]
  n_matched     <- sum(ft_sample_cols %in% valid_samples)
  match_pct     <- 100 * n_matched / length(valid_samples)

  message(sprintf("     → Sample matching: %d/%d metadata samples found in FT (%.1f%%)",
                  n_matched, length(valid_samples), match_pct))

  # If >=90% of metadata samples already match FT columns, no mapping needed
  if (match_pct >= 90) {
    return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                mapping_applied = FALSE))
  }

  # --- Detect SAMD (BioSample) vs Run accession mismatch ---
  # This is the most common ID mismatch in SRA-derived datasets: metadata uses
  # BioSample accessions (SAMD*) while feature tables use run accessions (SRR*)
  meta_has_samd       <- any(grepl("^SAMD", valid_samples))
  ft_has_run_accessions <- any(grepl("^(DRR|SRR|ERR)", ft_sample_cols))

  if (!meta_has_samd || !ft_has_run_accessions) {
    message("     ⚠️ Low match but no SAMD↔Run mismatch detected")
    return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                mapping_applied = FALSE))
  }

  message("     → Auto-detected SAMD↔Run ID mismatch. Attempting mapping...")

  # Reload the original (unfiltered) metadata to search for mapping columns
  orig_meta <- load_table_robust(orig_metadata, "SampleID")
  message(sprintf("       Metadata dimensions: %d rows × %d columns",
                  nrow(orig_meta), ncol(orig_meta)))

  # --- Strategy 1: Find mapping column by name ---
  # Look for columns named "run", "accession", "srr", etc.
  drr_col <- grep("run|accession|srr|drr|err|biosample.*run|experiment|library.*name",
                  names(orig_meta), value = TRUE, ignore.case = TRUE)[1]

  # --- Strategy 2: Find mapping column by content ---
  # Scan all columns for SRA run accession patterns
  if (is.na(drr_col)) {
    message("       → Searching columns for run accession codes...")
    for (col_name in names(orig_meta)[-1]) {
      values     <- as.character(orig_meta[[col_name]])
      n_run_codes <- sum(grepl("^(DRR|SRR|ERR)", values, ignore.case = TRUE), na.rm = TRUE)
      # Require at least 30% of samples to have run codes to avoid false positives
      if (n_run_codes > length(valid_samples) * 0.3) {
        drr_col <- col_name
        message(sprintf("       ✓ Found run codes in '%s' (%d matches)", col_name, n_run_codes))
        break
      }
    }
  } else {
    message(sprintf("       ✓ Using column '%s' for mapping", drr_col))
  }

  # --- Strategy 3: Sequential positional mapping (last resort) ---
  # Assumes samples are in the same order in metadata and feature table.
  # This is fragile and logs a prominent warning.
  if (is.na(drr_col)) {
    message("       ⚠️ No mapping column found. Attempting sequential mapping...")
    message("          WARNING: This assumes samples are in the same order!")

    if (nrow(meta_clean) <= length(ft_sample_cols)) {
      mapped_samples <- ft_sample_cols[1:nrow(meta_clean)]
      mapping_df     <- tibble(SAMD = meta_clean[[sample_id_col]], RUN = mapped_samples)
      message(sprintf("       ✓ Mapped %d samples sequentially", nrow(mapping_df)))
      message(sprintf("       Example: %s → %s", mapping_df$SAMD[1], mapping_df$RUN[1]))

      # Replace BioSample IDs with run accessions in metadata
      valid_samples <- mapping_df$RUN
      meta_clean <- meta_clean %>%
        mutate(original_id = .data[[sample_id_col]]) %>%
        left_join(mapping_df, by = c("original_id" = "SAMD")) %>%
        mutate(!!sample_id_col := coalesce(RUN, .data[[sample_id_col]])) %>%
        select(-RUN, -original_id)

      n_matched <- sum(ft_sample_cols %in% valid_samples)
      match_pct <- 100 * n_matched / length(valid_samples)
      message(sprintf("       → After sequential mapping: %d/%d matched (%.1f%%)",
                      n_matched, length(valid_samples), match_pct))
      return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                  mapping_applied = TRUE))
    } else {
      message("       ❌ Sequential mapping not possible (more meta samples than FT samples)")
      return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                  mapping_applied = FALSE))
    }
  }

  # --- Apply named-column mapping ---
  # Build a BioSample -> Run accession lookup table
  mapping_df <- orig_meta %>%
    select(SAMD = 1, RUN = all_of(drr_col)) %>%
    filter(!is.na(RUN), RUN != "", SAMD %in% valid_samples) %>%
    mutate(RUN = trimws(RUN))

  if (nrow(mapping_df) == 0) {
    message("       ❌ Mapping column found but no valid mappings extracted")
    return(list(valid_samples = valid_samples, meta_clean = meta_clean,
                mapping_applied = FALSE))
  }

  # Replace BioSample IDs with run accessions in metadata
  message(sprintf("       ✓ Created mapping for %d samples", nrow(mapping_df)))
  valid_samples <- mapping_df$RUN
  meta_clean <- meta_clean %>%
    left_join(mapping_df, by = setNames("SAMD", sample_id_col)) %>%
    mutate(!!sample_id_col := coalesce(RUN, .data[[sample_id_col]])) %>%
    select(-RUN)

  n_matched <- sum(ft_sample_cols %in% valid_samples)
  match_pct <- 100 * n_matched / length(valid_samples)
  message(sprintf("       → After column mapping: %d/%d matched (%.1f%%)",
                  n_matched, length(valid_samples), match_pct))

  return(list(valid_samples = valid_samples, meta_clean = meta_clean,
              mapping_applied = TRUE))
}

# ==============================================================================
# MAIN PROCESSING LOOP
# Iterates over each study in the manifest, applying the four-step cleaning
# pipeline: metadata filter -> taxonomy filter -> feature table sync -> FASTA sync
# ==============================================================================

for (i in seq_len(nrow(manifest))) {
  study      <- manifest[i, ]
  study_name <- study$study_name

  message(paste0("\n>>> Processing Study: ", study_name))

  # --- Parse Metadata Filtering Parameters ---
  # The manifest may optionally specify a metadata column and value(s) to filter
  # on (e.g., environment="marine,estuarine"). This enables subsetting studies
  # that span multiple environments to retain only the target habitat.
  has_metadata_var   <- "metadata_variable" %in% names(study)
  has_metadata_val   <- "metadata_value"    %in% names(study)
  should_filter_metadata <- FALSE
  TARGET_ENVIRONMENTS    <- NULL
  target_col             <- NULL

  if (has_metadata_var && has_metadata_val) {
    metadata_var_raw <- study$metadata_variable
    metadata_val_raw <- study$metadata_value
    var_is_valid <- !is.na(metadata_var_raw) && nchar(trimws(metadata_var_raw)) > 0
    val_is_valid <- !is.na(metadata_val_raw) && nchar(trimws(metadata_val_raw)) > 0
    if (var_is_valid && val_is_valid) {
      should_filter_metadata <- TRUE
      target_col          <- trimws(metadata_var_raw)
      # Comma-separated values allow multiple target environments
      TARGET_ENVIRONMENTS <- trimws(unlist(strsplit(as.character(metadata_val_raw), ",")))
      message(paste("    Metadata filtering: ENABLED"))
      message(paste("    Target variable:", target_col))
      message(paste("    Target values:", paste(TARGET_ENVIRONMENTS, collapse = ", ")))
    } else {
      message("    Metadata filtering: DISABLED (no variable/value specified)")
    }
  } else {
    message("    Metadata filtering: DISABLED (columns not in manifest)")
  }

  # --- File Existence Check ---
  # Verify all four input files exist before processing to fail early with
  # informative messages rather than cryptic downstream errors
  orig_counts   <- study$asv_counts_path
  orig_taxonomy <- study$taxonomy_path
  orig_metadata <- study$metadata_path
  orig_fasta    <- study$asv_fasta_path

  paths_orig <- c(counts = orig_counts, taxonomy = orig_taxonomy,
                  metadata = orig_metadata, fasta = orig_fasta)
  missing <- paths_orig[!file.exists(paths_orig)]

  if (length(missing) > 0) {
    message("  ⚠️ Skipping study due to missing files:")
    for (nm in names(missing)) message("     - ", nm, ": ", missing[[nm]])
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  # Derive output paths by inserting "clean/" subdirectory
  clean_counts   <- make_clean_path(orig_counts)
  clean_taxonomy <- make_clean_path(orig_taxonomy)
  clean_metadata <- make_clean_path(orig_metadata)
  clean_fasta    <- make_clean_path(orig_fasta)


  # Create output directories
for (p in unique(c(dirname(clean_counts), dirname(clean_taxonomy),
                   dirname(clean_metadata), dirname(clean_fasta)))) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

  # --- [1/4] Metadata Filtering ---
  # Load metadata and optionally subset to target environment(s)
  message("  [1/4] Processing Metadata...")
  meta_df       <- load_table_robust(orig_metadata, "SampleID")
  sample_id_col <- "SampleID"

  if (should_filter_metadata && !is.null(target_col) && target_col %in% names(meta_df)) {
    meta_clean <- meta_df %>% filter(.data[[target_col]] %in% TARGET_ENVIRONMENTS)
    message(sprintf("     - Filtered by %s: kept %d / %d samples",
                    target_col, nrow(meta_clean), nrow(meta_df)))
  } else {
    if (should_filter_metadata && !is.null(target_col) && !target_col %in% names(meta_df)) {
      message(sprintf("  ⚠️ Column '%s' not found in metadata — keeping all samples", target_col))
    }
    meta_clean <- meta_df
    message(paste("     - No filtering applied. Kept all", nrow(meta_clean), "samples"))
  }

  # Enforce minimum sample count to avoid underpowered studies
  if (nrow(meta_clean) < MIN_SAMPLES) {
    message(sprintf("  ⚠️ Only %d samples (minimum %d). Skipping.", nrow(meta_clean), MIN_SAMPLES))
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  # --- [2/4] Taxonomy Filtering ---
  # Remove Chloroplast and Mitochondria ASVs. These are eukaryotic organellar
  # 16S sequences that co-amplify with bacterial/archaeal targets and would
  # inflate diversity estimates if retained.
  message("  [2/4] Filtering Taxonomy...")
  tax_df    <- load_table_robust(orig_taxonomy, "ASV_ID")
  tax_col   <- names(tax_df)[2]
  tax_clean <- tax_df %>%
    filter(!grepl(TAXA_TO_REMOVE, .data[[tax_col]], ignore.case = TRUE))
  message(sprintf("     - Removed %d ASVs matching Chloroplast/Mitochondria",
                  nrow(tax_df) - nrow(tax_clean)))
  write_tsv(tax_clean, clean_taxonomy)
  valid_asvs <- tax_clean[["ASV_ID"]]

  # --- [3/4] Feature Table Cleaning ---
  # Sync feature table to: (a) retained ASVs from taxonomy filter, (b) retained
  # samples from metadata filter. Then remove zero-count ASVs and low-read samples.
  message("  [3/4] Cleaning Feature Table...")

  ft_df         <- load_feature_table(orig_counts)
  qiime_headers <- attr(ft_df, "qiime_headers")
  ft_sample_cols <- names(ft_df)[-1]

  message(sprintf("     → Feature table loaded: %d ASVs × %d samples",
                  nrow(ft_df), length(ft_sample_cols)))

  # Per-sample read totals BEFORE any filtering
  sample_totals_orig <- colSums(ft_df[, -1], na.rm = TRUE)
  zero_orig          <- sum(sample_totals_orig == 0)
  message(sprintf("     → Pre-filter read totals: min=%.0f, median=%.0f, max=%.0f",
                  min(sample_totals_orig), median(sample_totals_orig), max(sample_totals_orig)))
  if (zero_orig > 0) {
    message(sprintf("     ⚠️ %d/%d samples (%.1f%%) have ZERO reads before filtering",
                    zero_orig, length(sample_totals_orig),
                    100 * zero_orig / length(sample_totals_orig)))
  }

  # --- Reconcile sample IDs between feature table and metadata ---
  # Public datasets frequently have BioSample vs Run accession mismatches
  mapping_result <- auto_detect_and_map_samples(
    ft_sample_cols, meta_clean, sample_id_col, orig_metadata, study_name)
  valid_samples <- mapping_result$valid_samples
  meta_clean    <- mapping_result$meta_clean

  n_matched <- sum(ft_sample_cols %in% valid_samples)
  match_pct <- 100 * n_matched / length(valid_samples)

  # Require at least 50% sample overlap to proceed
  if (match_pct < 50) {
    message(sprintf("     ❌ Only %.1f%% of metadata samples matched in FT. Skipping.", match_pct))
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  # Subset to intersection of matched samples and valid (non-organellar) ASVs
  cols_to_keep <- c("ASV_ID", intersect(ft_sample_cols, valid_samples))
  if (length(cols_to_keep) == 1) {
    message("     ❌ FATAL: No matching sample columns after mapping. Skipping.")
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  ft_clean <- ft_df %>%
    select(all_of(cols_to_keep)) %>%
    filter(.data[["ASV_ID"]] %in% valid_asvs)

  # Remove ASVs with zero reads across retained samples
  # (these had all their reads in samples that were filtered out)
  asv_totals <- rowSums(ft_clean[, -1], na.rm = TRUE)
  ft_clean   <- ft_clean[asv_totals > 0, ]
  message(sprintf("     → After ASV zero-filter: %d ASVs retained", nrow(ft_clean)))

  # Per-sample totals AFTER ASV filter -- before sample removal
  sample_totals_post <- colSums(ft_clean[, -1], na.rm = TRUE)
  message(sprintf("     → Post-ASV-filter sample totals: min=%.0f, median=%.0f, max=%.0f",
                  min(sample_totals_post), median(sample_totals_post), max(sample_totals_post)))

  # Remove samples below read threshold
  empty_samples <- names(sample_totals_post)[sample_totals_post < MIN_READS_PER_SAMPLE]
  if (length(empty_samples) > 0) {
    message(sprintf("     ⚠️ Removing %d samples with < %d reads",
                    length(empty_samples), MIN_READS_PER_SAMPLE))
    message(sprintf("        Examples: %s", paste(head(empty_samples, 5), collapse = ", ")))
    # Show read totals for removed samples (diagnostic)
    removed_totals <- sort(sample_totals_post[empty_samples])
    message(sprintf("        Removed sample totals: %s",
                    paste(head(removed_totals, 10), collapse = ", ")))
    ft_clean <- ft_clean %>%
      select(all_of(c("ASV_ID", setdiff(names(ft_clean)[-1], empty_samples))))
  }

  message(sprintf("     → Final Feature Table: %d ASVs × %d samples",
                  nrow(ft_clean), ncol(ft_clean) - 1))

  # Synchronise metadata and FASTA to match the final feature table
  final_valid_asvs <- ft_clean[["ASV_ID"]]
  final_ft_samples <- names(ft_clean)[-1]

  meta_clean <- meta_clean %>% filter(.data[[sample_id_col]] %in% final_ft_samples)
  message(sprintf("     → Aligned metadata: %d samples", nrow(meta_clean)))

  # Write metadata
  write_tsv(meta_clean, clean_metadata)

  # Write feature table (preserve QIIME2 header if present)
  # QIIME2 format requires the original comment line and a #-prefixed header row
  if (!is.null(qiime_headers)) {
    writeLines(c(qiime_headers[1],
                 paste0("#", paste(colnames(ft_clean), collapse = "\t"))),
               clean_counts)
    write_tsv(ft_clean, clean_counts, append = TRUE, col_names = FALSE)
  } else {
    write_tsv(ft_clean, clean_counts, col_names = TRUE)
  }

  # --- [4/4] FASTA Synchronisation ---
  # Subset representative sequences to match the cleaned feature table ASVs
  message("  [4/4] Cleaning FASTA...")
  seqs       <- readDNAStringSet(orig_fasta)
  seqs_clean <- seqs[names(seqs) %in% final_valid_asvs]
  message(sprintf("     - Removed %d sequences (%d retained)",
                  length(seqs) - length(seqs_clean), length(seqs_clean)))
  writeXStringSet(seqs_clean, clean_fasta)

  # --- Final Validation ---
  # Guard against studies that pass all filters but end up empty
  if (length(seqs_clean) == 0 || nrow(ft_clean) == 0 || ncol(ft_clean) <= 1) {
    message("  ⚠️ No valid ASVs or samples remaining. Excluding from manifest.")
    manifest_clean[i, c("asv_counts_path","taxonomy_path","metadata_path","asv_fasta_path")] <- NA_character_
    next
  }

  # Update manifest to point to cleaned files
  manifest_clean$asv_counts_path[i] <- clean_counts
  manifest_clean$taxonomy_path[i]   <- clean_taxonomy
  manifest_clean$metadata_path[i]   <- clean_metadata
  manifest_clean$asv_fasta_path[i]  <- clean_fasta

  message("  ✅ Study cleaned and validated")
}

# ==============================================================================
# POST-CLEANING VALIDATION
# Re-load each cleaned study and verify that feature table sample IDs perfectly
# match metadata sample IDs. This catches any edge cases where the cleaning
# logic produced misaligned outputs.
# ==============================================================================

message("\n========================================================")
message("POST-CLEANING VALIDATION")
message("========================================================")

# Remove studies that failed cleaning (NA paths)
manifest_clean_final <- manifest_clean %>%
  filter(!is.na(asv_counts_path) & !is.na(asv_fasta_path))

validation_issues <- 0

for (i in seq_len(nrow(manifest_clean_final))) {
  study      <- manifest_clean_final[i, ]
  study_name <- study$study_name

  clean_ft   <- load_feature_table(study$asv_counts_path)
  clean_meta <- load_table_robust(study$metadata_path, "SampleID")

  ft_samples   <- names(clean_ft)[-1]
  meta_samples <- clean_meta$SampleID
  n_match      <- sum(ft_samples %in% meta_samples)
  match_pct    <- 100 * n_match / length(ft_samples)

  if (match_pct < 100) {
    message(sprintf("  %s: ⚠️ %d/%d FT samples in metadata (%.1f%%)",
                    study_name, n_match, length(ft_samples), match_pct))
    missing_from_meta <- setdiff(ft_samples, meta_samples)
    extra_in_meta     <- setdiff(meta_samples, ft_samples)
    if (length(missing_from_meta) > 0)
      message(sprintf("    Missing from metadata: %s", paste(head(missing_from_meta, 5), collapse = ", ")))
    if (length(extra_in_meta) > 0)
      message(sprintf("    Extra in metadata: %s",    paste(head(extra_in_meta, 5), collapse = ", ")))
    validation_issues <- validation_issues + 1
  } else {
    message(sprintf("  %s: ✅ Perfect match (%d samples)", study_name, n_match))
  }
}

# --- Write Final Manifest ---
write_tsv(manifest_clean_final, MANIFEST_OUT)

message("\n========================================================")
message(paste("Cleaning complete. Updated manifest:", MANIFEST_OUT))
message(paste("Studies in final manifest:", nrow(manifest_clean_final)))
if (validation_issues > 0) {
  message(sprintf("⚠️ %d studies have validation issues — review above", validation_issues))
} else {
  message("✅ All studies validated successfully!")
}
message("Use secat_manifest_clean.tsv for downstream MESAP runs.")
message("========================================================")
