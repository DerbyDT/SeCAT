#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   06_analyse_real.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 6 — Real Data Trimming Analysis (Worker)
# PURPOSE:  Progressively trim one study's sequences and measure beta-diversity
#           degradation at each step to build the real degradation curve.
#
# OVERVIEW:
#   This is the empirical counterpart to the simulation worker (05_sim_worker.R).
#   For a single study, it:
#     1. Loads the study's aligned ASV sequences and feature table
#     2. Determines the number of trim steps needed to reach (and pass) the
#        consensus region boundary
#     3. At each step, trims sequences from the amplicon edges inward, re-
#        clusters with VSEARCH, and computes Bray-Curtis dissimilarity relative
#        to the untrimmed baseline
#     4. Identifies changepoints (PELT) in the degradation curve
#     5. Tracks taxonomic retention and core taxa dynamics
#
#   The resulting degradation curves are compared against simulation null
#   distributions in 07_aggregate.R to determine whether the study can be
#   safely trimmed to the consensus region.
#
# INPUTS:
#   - Command-line argument: study_name (string)
#   - SECAT_MANIFEST_PATH: cleaned manifest TSV
#   - output/intermediate/study_alignment_coords.csv (study coordinates)
#   - output/intermediate/aligned_fastas/<study>_aligned.fasta
#   - output/intermediate/consensusregioninfo.csv (consensus boundaries)
#
# OUTPUTS:
#   - output/real_data_results/<study>/<study>_results.rds
#     (list containing: dissimilarity curves, retention data, changepoints,
#      OTU tables per taxonomic level, taxon impact analysis, trim metadata)
#
# DEPENDENCIES:
#   - R/secat_config.R: global paths and configuration
#   - R/secat_utils.R: run_trim_analysis(), calculate_dissimilarity(),
#     calculate_taxonomic_retention(), calculate_changepoint_thresholds(),
#     analyze_core_taxa(), analyze_taxon_impact(), sync_data_for_study()
#   - R/secat_consensus.R: find_largest_overlapping_clique()
#
# CALLED BY:
#   - modules/local/analyse_real.nf (ANALYSE_REAL process)
# ==============================================================================

#===============================================================================
# SECTION 1: SETUP
#===============================================================================

cat("--- Loading Libraries & Helper Functions ---\n")
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
})

source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R"))

# Standardize step mode (Global setting)
current_step_mode <- if(exists("TRIM_STEP_MODE")) TRIM_STEP_MODE else "absolute"

# Validate configuration
if (!current_step_mode %in% c("scaled", "absolute")) {
    warning(sprintf("Invalid TRIM_STEP_MODE '%s'. Defaulting to 'absolute'.", current_step_mode))
    current_step_mode <- "absolute"
}

# Load core engines
source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_utils.R"))

# --- Get Task ID / Study Name ---
# Fetch arguments explicitly passed from the Nextflow command line
args <- commandArgs(trailingOnly = TRUE)

# Safety check: ensure Nextflow actually passed an argument
if (length(args) < 1) {
  stop("FATAL: Requires 'study_name' as a command-line argument.")
}

# The first argument passed will be our target study
target_study <- args[1]

cat("--- Loading Data ---\n")
# Load manifest to identify the study
manifest <- readr::read_tsv(SECAT_MANIFEST_PATH, show_col_types = FALSE)

# Filter the manifest to find the row matching our target study
job_info <- subset(manifest, study_name == target_study)

# Safety check: ensure we found exactly one matching study
if (nrow(job_info) != 1) {
  stop(paste("FATAL: Could not find exactly one entry for study:", target_study))
}

# Extract the variables exactly as before
study_name <- job_info$study_name
primer_name_from_manifest <- job_info$primer_name
cat(paste("Processing study:", study_name, "| Primer set:", primer_name_from_manifest, "\n"))

#===============================================================================
# REPRODUCIBILITY: Study-based deterministic seeding (FIXED)
#===============================================================================

if (!requireNamespace("digest", quietly = TRUE)) {
  message("  -> Installing 'digest' package...")
  install.packages("digest", repos = "https://cloud.r-project.org")
}

# Generate seed safely (always fits in 32-bit signed int)
study_hash <- digest::digest(study_name, algo = "md5")

# Take 6 hex chars instead of 8 (max value: 16,777,215)
study_seed <- strtoi(substr(study_hash, 1, 6), 16L)


set.seed(study_seed)
cat(sprintf("  -> [REPRODUCIBILITY] Seed: %d (Study: %s)\n", study_seed, study_name))

# ===================================================================
# SECTION 2: COORDINATE & CONSENSUS LOGIC
# ===================================================================
cat("--- Looking up coordinates and calculating consensus region ---\n")

coords_for_run <- NULL
consensus_start_global <- NULL
consensus_end_global <- NULL
USE_ALIGNMENT_LOGIC <- (ANALYSIS_MODE == "study")

# Load the clique consensus function
source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R", "secat_consensus.R"))

if (ANALYSIS_MODE == "primer") {
    # ---------------------------------------------------------
    # PRIMER MODE: Theoretical Clique Consensus
    # ---------------------------------------------------------
    coords_path <- file.path("output/intermediate/primer_coords_phase1_output.csv")
    if (!file.exists(coords_path)) stop("[FATAL] primer_coords_phase1_output.csv not found.")

    # 1. Load Primer Database (Dictionary of Primer Name -> Coords)
    primer_db <- readr::read_csv(coords_path, show_col_types = FALSE)

    # Normalize columns if needed (handle your underscore naming issue here too)
    if ("primerstart" %in% names(primer_db)) {
        primer_db <- primer_db %>% dplyr::rename(primer_start = primerstart, primer_end = primerend)
    }

    # 2. Get specific coords for THIS study (to calculate distance later)
    coords_for_run <- primer_db %>% dplyr::filter(primer_name == primer_name_from_manifest)

    if (nrow(coords_for_run) == 0) {
        stop(paste0("[FATAL] Primer '", primer_name_from_manifest, "' not found in primer database."))
    }

    # 3. CALCULATE CLIQUE CONSENSUS (The Upgrade)
    # We join the manifest with the primer DB to get theoretical coords for ALL studies
    cat("   -> Calculating theoretical consensus across all studies in manifest...\n")

    # Ensure manifest has clean primer names matching the DB
    manifest_theoretical <- manifest %>%
        dplyr::inner_join(primer_db, by = c("primer_name" = "primer_name"))

    if (nrow(manifest_theoretical) == 0) {
        # Fallback if join fails entirely (e.g. names don't match)
        cat("   -> [WARN] Manifest primer names don't match DB. Using single primer consensus.\n")
        consensus_start_global <- coords_for_run$primer_start[1]
        consensus_end_global <- coords_for_run$primer_end[1]
    } else {
        # Run the clique algorithm on the theoretical coordinates
        consensus_result <- find_largest_overlapping_clique(
            starts = manifest_theoretical$primer_start,
            ends = manifest_theoretical$primer_end,
            study_names = manifest_theoretical$study_name,
            min_overlap = 50
        )

        consensus_start_global <- consensus_result$start
        consensus_end_global <- consensus_result$end
    }

} else if (ANALYSIS_MODE == "study") {
    # ---------------------------------------------------------
    # STUDY MODE: Empirical Clique Consensus (Existing Logic)
    # ---------------------------------------------------------
    cat("   -> In study mode. Calculating global consensus from empirical alignments.\n")
    coords_path <- file.path("output/intermediate/study_alignment_coords.csv")

    if (!file.exists(coords_path)) stop("[FATAL] study_alignment_coords.csv not found.")

    all_coords_data <- readr::read_csv(coords_path, show_col_types = FALSE)

    # Run clique algorithm on empirical coordinates
    consensus_result <- find_largest_overlapping_clique(
        starts = all_coords_data$ref_start,
        ends = all_coords_data$ref_end,
        study_names = all_coords_data$study_name,
        min_overlap = 50
    )

    consensus_start_global <- consensus_result$start
    consensus_end_global <- consensus_result$end

    # Get specific coords for THIS study
    coords_for_run <- all_coords_data %>% dplyr::filter(study_name == !!study_name)

    if (nrow(coords_for_run) != 1) stop("[FATAL] Coordinate lookup failed for this study.")
}

cat(sprintf("   -> Consensus Region: %d - %d\n", consensus_start_global, consensus_end_global))

# ---------------------------------------------------------
# Define Study Boundaries (Standardized)
# ---------------------------------------------------------
# Now we can safely extract start/end regardless of mode
study_s <- if (ANALYSIS_MODE == "primer") coords_for_run$primer_start[1] else coords_for_run$ref_start[1]
study_e <- if (ANALYSIS_MODE == "primer") coords_for_run$primer_end[1] else coords_for_run$ref_end[1]

# ===================================================================
# SECTION 3: DATA PREPARATION
# ===================================================================

# --- Load Original Study Data ---
cat("Step 1: Synchronizing study data...\n")
study_output_dir <- file.path(REAL_DATA_RESULTS_DIR, study_name)
dir.create(study_output_dir, recursive = TRUE, showWarnings = FALSE)

# sync_data_for_study loads the OTU table and Sequences, checking for ID mismatches
study_data <- sync_data_for_study(job_info)
otu_table <- study_data$otu_table
sequences <- study_data$sequences
has_taxonomy <- study_data$has_taxonomy

# --- Load Pre-Aligned Sequences (if applicable) ---
cat(paste("Step 2: Preparing trim analysis (Use Alignment =", USE_ALIGNMENT_LOGIC, ")\n"))

aligned_seqs_for_run <- NULL

if (USE_ALIGNMENT_LOGIC) {
    # Alignment object generated in Phase 2
    align_path <- file.path("output", "intermediate", "aligned_fastas", paste0(study_name, "_aligned.fasta"))
    if (file.exists(align_path)) {
        cat("  -> Loading alignment object from Step 02...\n")
        aligned_seqs_for_run <- Biostrings::readDNAStringSet(align_path)
    } else {
        cat("  -> WARNING: Alignment file not found. Falling back to physical trimming.\n")
    }
}

# ===================================================================
# SECTION 4: TRIM PARAMETER CALCULATION
# ===================================================================

# --- Define Boundaries ---
study_s <- if (ANALYSIS_MODE == "primer") coords_for_run$primer_start else coords_for_run$ref_start
study_e <- if (ANALYSIS_MODE == "primer") coords_for_run$primer_end else coords_for_run$ref_end

# --- Check Outlier Status ---
safe_consensus_start <- consensus_start_global
safe_consensus_end <- consensus_end_global
is_outlier <- FALSE

if (safe_consensus_end <= safe_consensus_start || study_s >= safe_consensus_end || study_e <= safe_consensus_start) {
    cat("  -> WARNING: Study is an OUTLIER (does not overlap consensus).\n")
    safe_consensus_start <- study_s
    safe_consensus_end <- study_e
    is_outlier <- TRUE
}

# --- Calculate Steps & Increment ---
# Matches Phase 3 Logic Exactly
final_num_steps <- 0
final_increment <- 0

# Calculate Total Distance to Consensus
dist_start <- 0
dist_end <- 0
if (study_s < safe_consensus_start) dist_start <- safe_consensus_start - study_s
if (study_e > safe_consensus_end)   dist_end <- study_e - safe_consensus_end
total_dist_cols <- dist_start + dist_end

buffer_steps <- if(exists("CONSENSUS_BUFFER_STEPS")) CONSENSUS_BUFFER_STEPS else 10
current_step_mode <- if(exists("TRIM_STEP_MODE")) TRIM_STEP_MODE else "absolute"

if (current_step_mode == "absolute") {
    # --- ABSOLUTE RESOLUTION MODE (Fixed Increment) ---
    target_increment <- if (exists("TRIM_INCREMENT")) TRIM_INCREMENT else 5

    if (total_dist_cols > 0) {
        steps_to_consensus <- ceiling(total_dist_cols / target_increment)
        cat(sprintf("  -> Mode: ABSOLUTE. Dist: %d. Inc: %d.\n", total_dist_cols, target_increment))
        cat(sprintf("  -> Steps required: %d\n", steps_to_consensus))

        final_num_steps <- steps_to_consensus + buffer_steps
        final_increment <- target_increment
    } else {
        final_num_steps <- buffer_steps
        final_increment <- target_increment
    }

    # Safety Cap
    safety_limit <- if(exists("MAX_ABSOLUTE_TRIM_STEPS")) MAX_ABSOLUTE_TRIM_STEPS else 2000
    if (final_num_steps > safety_limit) {
        final_num_steps <- safety_limit
    }

} else {
    # --- SCALED MODE (Fixed Steps) ---
    target_steps <- if (exists("DEFAULT_MAX_TRIM_STEPS")) DEFAULT_MAX_TRIM_STEPS else 40

    if (total_dist_cols > 0) {
        dynamic_inc <- ceiling(total_dist_cols / target_steps)
        cat(sprintf("  -> Mode: SCALED. Dist: %d. Inc: %d.\n", total_dist_cols, dynamic_inc))

        final_increment <- dynamic_inc
        final_num_steps <- target_steps + buffer_steps
    } else {
        final_increment <- if (exists("TRIM_INCREMENT")) TRIM_INCREMENT else 5
        final_num_steps <- buffer_steps
    }
}

# ===================================================================
# SECTION 5: EXECUTION (TRIMMING & METRICS)
# ===================================================================

# --- Execute Trimming Engine ---
cat("Step 2: Executing trim analysis...\n")
trim_results <- run_trim_analysis(
    sequences = sequences,
    vsearch_path = VSEARCH_PATH,
    output_dir = study_output_dir,
    num_steps = final_num_steps,
    mode = "both",
    increment = final_increment,
    primer_start    = study_s,
    primer_end      = study_e,
    consensus_start = safe_consensus_start,
    consensus_end   = safe_consensus_end,
    use_alignment   = USE_ALIGNMENT_LOGIC,
    aligned_sequences = aligned_seqs_for_run
)
trim_results$is_outlier <- is_outlier

# --- Analyze Results ---
cat("Step 3: Calculating metrics...\n")

# 1. Build OTU Tables
otu_tables_per_level <- analyze_all_taxonomic_levels(otu_table, trim_results$uc_files)

# 2. Max Valid Step
max_valid_trim_step <- calculate_max_valid_trim(otu_tables_per_level, increment = final_increment)

# 3. Dissimilarity (Beta Diversity)
cat("Step 4: Calculating dissimilarity...\n")
dissim_data <- calculate_dissimilarity(
  otu_tables_per_level,
  final_num_steps,
  increment = final_increment
)

# 4. Retention (Alpha-like)
retention_data <- calculate_taxonomic_retention(
  otu_tables_per_level,
  final_num_steps,
  increment = final_increment
)

# 5. Changepoints (Tipping Points)
cat("   -> Calculating changepoint thresholds...\n")
threshold_data <- calculate_changepoint_thresholds(dissim_data)

# 6. Core Taxa & Impacts (Biological Insight)
cat("   -> Analyzing core taxa dynamics...\n")
use_prevalence <- if(exists("CORE_PREVALENCE_THRESHOLD")) CORE_PREVALENCE_THRESHOLD else 0.75
use_top_n <- if(exists("TOP_N_CORE_TAXA")) TOP_N_CORE_TAXA else 10

core_taxa_data <- analyze_core_taxa(
  otu_tables_per_level,
  prevalence_threshold = use_prevalence,
  top_n = use_top_n
)

# Tracks every taxon's trajectory for plots
impacted_taxa_data <- analyze_taxon_impact(
  otu_tables_per_level,
  num_steps = final_num_steps,
  increment = final_increment
)

taxon_impacts_for_plots <- list(
  core_taxa = core_taxa_data,
  impacted_taxa = impacted_taxa_data$impacted_taxa
)

if (is.null(core_taxa_data)) {
  cat("   ! Warning: No core taxa identified.\n")
} else {
  cat(sprintf("   + Identified dynamics for %d core taxa.\n", nrow(core_taxa_data)))
}

# ===================================================================
# SECTION 6: SAVE OUTPUT
# ===================================================================
cat("Step 5: Saving results...\n")

if (!exists("output_dir")) {
  output_dir <- file.path("output", "real_data_results", study_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

final_output <- list(
  study_name = study_name,
  primer_name = primer_name_from_manifest,
  has_taxonomy = has_taxonomy,
  target_trim_bp = trim_results$total_trim_to_consensus,
  head_proportion = trim_results$head_proportion,
  num_steps_used = final_num_steps,
  increment = final_increment,
  max_valid_trim_step = max_valid_trim_step,
  otu_table = otu_table,
  otu_tables_per_level = otu_tables_per_level,
  uc_files = trim_results$uc_files,
  dissim_data = dissim_data,
  retention_data = retention_data,
  threshold_data = threshold_data,
  taxon_impacts = taxon_impacts_for_plots,
  analysis_mode = ANALYSIS_MODE,
  use_alignment = USE_ALIGNMENT_LOGIC,
  trim_step_mode = current_step_mode,
  is_outlier = is_outlier
)

save_path <- file.path(output_dir, paste0(study_name, "_results.rds"))
saveRDS(final_output, save_path)
cat(sprintf("Saved results to: %s\n", save_path))

cat("--- Worker Script Completed Successfully ---\n")
