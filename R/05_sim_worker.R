#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   05_sim_worker.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# STAGE:    Stage 5 - Simulation Worker (parallel execution)
# PURPOSE:  Generate one replicate of the null-model degradation curve
#
# OVERVIEW:
#   Each invocation of this script produces a single simulation replicate for
#   one task (study or primer). It draws a synthetic microbial community from
#   the SILVA reference subset, applies realistic PCR amplification bias
#   (GC-content dependent), Illumina-style position-dependent sequencing
#   errors, and optional chimera formation (inspired by Grinder; Angly et al.
#   2012, NAR 40:e94). The simulated amplicons are then progressively trimmed
#   from the consensus region boundaries inward, reclustered with VSEARCH at
#   each step, and Bray-Curtis dissimilarity is measured against the untrimmed
#   baseline. This builds the null distribution: "how much beta-diversity
#   change is expected from trimming alone, absent real biological signal?"
#
# INPUTS:
#   - output/intermediate/simulation_tasks.csv        Task manifest from Stage 4
#   - output/intermediate/simulation_reference_subset.fasta  SILVA subset
#   - output/intermediate/study_alignment_coords.csv  (study mode)
#   - output/intermediate/primer_coords_phase1_output.csv (primer mode)
#   - output/intermediate/consensusregioninfo.csv     Consensus region bounds
#   - Command-line args: TASK_ID, SEED
#
# OUTPUTS:
#   - output/simulation_results/{task_id}/seed_{seed}/results.rds
#     Contains: dissimilarity curves, retention metrics, changepoints,
#     taxon impact data, and simulation parameters for reproducibility
#
# DEPENDENCIES:
#   - R packages: dplyr, tidyr, ggplot2, readr, stringr, purrr, tibble,
#                 here, Biostrings
#   - Sourced: R/secat_config.R, R/secat_utils.R (trimming engine, community
#              generator, dissimilarity calculator, VSEARCH wrapper)
#   - External: VSEARCH binary (path from config)
#
# CALLED BY:
#   - Nextflow process: run_simulation (scattered across task_id x seed)
# ==============================================================================

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

#===============================================================================
# SECTION 1: SETUP & CONFIGURATION
#===============================================================================

# --- Load pipeline-wide configuration and utility functions ---
source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R"))
source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_utils.R"))

# Determine trim step mode: "scaled" uses a fixed number of steps with variable
# increment; "absolute" uses a fixed bp increment with variable step count.
current_step_mode <- if(exists("TRIM_STEP_MODE")) TRIM_STEP_MODE else "absolute"

if (!current_step_mode %in% c("scaled", "absolute")) {
    warning(sprintf("Invalid TRIM_STEP_MODE '%s'. Defaulting to 'absolute'.", current_step_mode))
    current_step_mode <- "absolute"
}

# --- Parse command-line arguments ---
# Nextflow passes TASK_ID (study/primer name) and SEED (replicate number)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("FATAL: Requires TASK_ID and SEED")

TASK_ID <- args[1]
SEED <- as.integer(args[2])
set.seed(SEED)

message(paste("Mode:", ANALYSIS_MODE, "| Task:", TASK_ID, "| Seed:", SEED))

# --- Load this worker's task specification from the master task list ---
tasks_path <- file.path(OUTDIR, "intermediate/simulation_tasks.csv")
if (!file.exists(tasks_path)) stop("FATAL: simulation_tasks.csv missing.")
all_tasks <- readr::read_csv(tasks_path, show_col_types = FALSE)
task_info <- all_tasks %>% dplyr::filter(task_id == TASK_ID, simulation_seed == SEED)
if (nrow(task_info) == 0) stop("FATAL: Task info missing.")

#===============================================================================
# SECTION 2: LOAD SILVA REFERENCE DATABASE (USE PREBUILT SUBSET)
#===============================================================================
# The full SILVA database is too large for each worker to load independently.
# Instead, 04_prepare_sims.R created a shared random subset that all workers
# draw from. This subset has already been cleaned (U->T conversion, gap
# standardisation) and is stored as a standard DNA FASTA.

message("  -> Loading SILVA reference database for null model simulation...")

subset_path <- file.path(OUTDIR, "intermediate/simulation_reference_subset.fasta")

if (!file.exists(subset_path)) {
  stop("FATAL: Simulation reference subset not found.\n",
       "Expected: ", subset_path, "\n",
       "This file should be created by scripts/04_prepare_sims.R\n",
       "Run: Rscript scripts/04_prepare_sims.R")
}

file_size_mb <- file.info(subset_path)$size / 1024^2
message(sprintf("  -> Using pre-built subset: %s (%.1f MB)",
                basename(subset_path), file_size_mb))

# Load the pre-cleaned subset (RNA->DNA conversion already applied)
ref_db <- tryCatch({
  Biostrings::readDNAStringSet(subset_path, format = "fasta", use.names = TRUE)
}, error = function(e) {
  stop("FATAL: Could not load simulation reference subset.\n",
       "Path: ", subset_path, "\n",
       "Error: ", conditionMessage(e), "\n",
       "Try regenerating with: Rscript scripts/04_prepare_sims.R")
})

# Sanity checks on loaded reference
n_seqs <- length(ref_db)
avg_width <- mean(width(ref_db))
message(sprintf("  -> Loaded %d sequences successfully", n_seqs))
message(sprintf("  -> Average sequence length: %.0f bp", avg_width))

if (n_seqs < 50) {
  warning("Very few sequences in simulation subset (", n_seqs, "). ",
          "Consider increasing SIMULATION_MAX_SILVA_SUBSET in config.")
}

#==============================================================================
# SECTION 3: LOAD COORDINATES (MODE-AWARE)
#==============================================================================
# Load the amplicon coordinates for this task and the global consensus region.
# These define where trimming starts and how many bases extend beyond the
# consensus (the "non-shared" flanking region that gets trimmed first).

consensus_start_global <- NULL
consensus_end_global <- NULL

if (ANALYSIS_MODE == "study") {

  # --- STUDY MODE: empirical coordinates from SILVA alignment ---
  coords_path <- file.path(OUTDIR, "intermediate/study_alignment_coords.csv")

  if (!file.exists(coords_path)) {
    stop("FATAL: Study coordinates file missing.\n",
         "Expected: ", coords_path, "\n",
         "In study mode, this file is created by scripts/02_study_mapping.R\n",
         "Run: qsub -t 1-N scripts/02_study_mapping.sge")
  }

  all_coords <- readr::read_csv(coords_path, show_col_types = FALSE)
  task_coords <- all_coords %>% dplyr::filter(study_name == TASK_ID)

  if (nrow(task_coords) == 0) {
    stop("FATAL: No coordinates found for study: ", TASK_ID)
  }

  study_s <- task_coords$ref_start
  study_e <- task_coords$ref_end

  message(sprintf("  -> Study Mode: Loaded coordinates for %s", TASK_ID))

  # Load global consensus region (computed by 04_prepare_sims.R)
  consensus_path <- file.path(OUTDIR, "intermediate/consensusregioninfo.csv")

  if (file.exists(consensus_path)) {
    consensus_info <- readr::read_csv(consensus_path, show_col_types = FALSE)
    consensus_start_global <- consensus_info$ConsensusStart[1]
    consensus_end_global <- consensus_info$ConsensusEnd[1]
  } else {
    stop("FATAL: consensusregioninfo.csv not found.\n",
         "This file is created by scripts/04_prepare_sims.R\n",
         "Run: Rscript scripts/04_prepare_sims.R")
  }

} else if (ANALYSIS_MODE == "primer") {

  # --- PRIMER MODE: theoretical coordinates from primer binding sites ---
  coords_path <- file.path(OUTDIR, "intermediate/primer_coords_phase1_output.csv")

  if (!file.exists(coords_path)) {
    stop("FATAL: Primer coordinates file missing.\n",
         "Expected: ", coords_path, "\n",
         "In primer mode, this file is created by scripts/01_primer_mapping.R\n",
         "Run: Rscript scripts/01_primer_mapping.R")
  }

  all_coords <- readr::read_csv(coords_path, show_col_types = FALSE)
  task_coords <- all_coords %>% dplyr::filter(primer_name == TASK_ID)

  if (nrow(task_coords) == 0) {
    stop("FATAL: No coordinates found for primer: ", TASK_ID)
  }

  # Primer mode uses primer binding positions as amplicon boundaries
  study_s <- task_coords$primer_start
  study_e <- task_coords$primer_end

  # Consensus already computed during Phase 1 and stored in primer coords
  consensus_start_global <- task_coords$consensus_start
  consensus_end_global <- task_coords$consensus_end

  message(sprintf("  -> Primer Mode: Loaded coordinates for %s", TASK_ID))

} else {
  stop("FATAL: Invalid ANALYSIS_MODE: ", ANALYSIS_MODE,
       ". Must be 'study' or 'primer'.")
}

# --- Validate consensus coordinates ---
if (is.null(consensus_start_global) || is.na(consensus_start_global) ||
    is.null(consensus_end_global) || is.na(consensus_end_global)) {
  stop("FATAL: Failed to load consensus region coordinates for task: ", TASK_ID)
}

message(sprintf("  -> Global Consensus Region: %d - %d",
               consensus_start_global, consensus_end_global))
message(sprintf("  -> Task Boundaries: %d - %d", study_s, study_e))

# Simulations always use alignment-based (columnar) trimming on SILVA coordinates
USE_ALIGNMENT <- TRUE

#===============================================================================
# SECTION 4: COMMUNITY GENERATION (GRINDER-INSPIRED)
#===============================================================================
# Generate a synthetic microbial community with realistic properties:
#   - Species drawn randomly from the SILVA subset
#   - Relative abundances follow a log-normal distribution (typical of real communities)
#   - PCR amplification bias modelled as GC-content dependent efficiency
#   - Illumina sequencing errors with position-dependent error profiles
#   - Optional chimera formation at a configurable rate
# This synthetic community serves as the "known truth" against which trimming
# effects are measured. Because we control the community composition, any
# observed beta-diversity change must be attributable to trimming artefacts.

message("Step 1: Simulating realistic community from SILVA...")

n_taxa <- if(exists("SIMULATION_COMMUNITY_SIZE")) SIMULATION_COMMUNITY_SIZE else 100

sim_community <- get_community(
  db_seq = ref_db,
  ntaxas = n_taxa,
  seed = SEED
)

sim_sequences_raw <- sim_community$sequences
otu_table <- sim_community$table

# The raw sequences retain SILVA alignment gaps for coordinate-based trimming
aligned_sim_seqs <- sim_sequences_raw  # Already aligned (from SILVA)

# VSEARCH requires ungapped sequences for de novo clustering
sim_sequences <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(sim_sequences_raw)))
names(sim_sequences) <- names(sim_sequences_raw)

message(sprintf("     -> Sequences in global coordinate system (Width: %d)",
                width(aligned_sim_seqs)[1]))

#===============================================================================
# SECTION 5: TRIM PARAMETER CALCULATION (MIRRORS REAL DATA)
#===============================================================================
# The trim parameters must exactly mirror those used in 06_analyse_real.R so
# the null distribution and real degradation curves are directly comparable.
# Key variables:
#   - dist_start/dist_end: how many alignment columns this study's amplicon
#     extends beyond the consensus region on each side
#   - total_dist_cols: total non-shared flanking region (trimmed first)
#   - buffer_steps: additional steps that trim INTO the consensus region to
#     characterise sensitivity beyond the consensus boundary

safe_consensus_start <- consensus_start_global
safe_consensus_end <- consensus_end_global
is_outlier <- FALSE

# Detect outlier studies whose amplicons do not overlap the consensus at all
if (safe_consensus_end <= safe_consensus_start ||
    study_s >= safe_consensus_end ||
    study_e <= safe_consensus_start) {

    message("  -> WARNING: Study is an OUTLIER (doesn't overlap global consensus).")
    safe_consensus_start <- study_s
    safe_consensus_end <- study_e
    is_outlier <- TRUE
}

# --- Calculate distance (in alignment columns) from study edges to consensus ---
dist_start <- 0
dist_end <- 0
if (study_s < safe_consensus_start) dist_start <- safe_consensus_start - study_s
if (study_e > safe_consensus_end)   dist_end <- study_e - safe_consensus_end
total_dist_cols <- dist_start + dist_end

buffer_steps <- if(exists("CONSENSUS_BUFFER_STEPS")) CONSENSUS_BUFFER_STEPS else 10

final_num_steps <- 0
final_increment <- 0

if (current_step_mode == "absolute") {
    # --- ABSOLUTE MODE: fixed bp increment, variable number of steps ---
    target_increment <- if (exists("TRIM_INCREMENT")) TRIM_INCREMENT else 5

    if (total_dist_cols > 0) {
        steps_to_consensus <- ceiling(total_dist_cols / target_increment)
        message(sprintf("  -> Mode: ABSOLUTE. Dist: %d. Inc: %d.", total_dist_cols, target_increment))

        # Steps to reach consensus + buffer steps beyond it
        final_num_steps <- steps_to_consensus + buffer_steps
        final_increment <- target_increment
    } else {
        # Study already within consensus; only buffer steps needed
        final_num_steps <- buffer_steps
        final_increment <- target_increment
    }

    # Safety cap to prevent runaway computation
    safety_limit <- if(exists("MAX_ABSOLUTE_TRIM_STEPS")) MAX_ABSOLUTE_TRIM_STEPS else 2000
    if (final_num_steps > safety_limit) {
        final_num_steps <- safety_limit
    }

} else {
    # --- SCALED MODE: fixed number of steps, variable bp increment ---
    # The increment is dynamically sized so that the target number of steps
    # spans exactly the non-shared flanking region.
    target_steps <- if (exists("DEFAULT_MAX_TRIM_STEPS")) DEFAULT_MAX_TRIM_STEPS else 40

    if (total_dist_cols > 0) {
        dynamic_inc <- ceiling(total_dist_cols / target_steps)
        final_increment <- dynamic_inc
        final_num_steps <- target_steps + buffer_steps
    } else {
        final_increment <- if (exists("TRIM_INCREMENT")) TRIM_INCREMENT else 5
        final_num_steps <- buffer_steps
    }
}

#===============================================================================
# SECTION 6: EXECUTION (TRIMMING & METRICS)
#===============================================================================
# Progressive trimming loop: at each step, remove `increment` more bases from
# both ends of the amplicon (outside the consensus inward), recluster with
# VSEARCH, and measure how much the resulting community profile has changed
# relative to the untrimmed baseline (step 0).

message(sprintf("Step 2: Running trim... Steps: %d, Inc: %d", final_num_steps, final_increment))

output_dir <- file.path(SIMULATION_RESULTS_DIR, TASK_ID, paste0("seed_", SEED))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Run the core trim-and-cluster engine (shared with 06_analyse_real.R)
trim_results <- run_trim_analysis(
    sequences = sim_sequences,
    vsearch_path = VSEARCH_PATH,
    output_dir = output_dir,
    num_steps = final_num_steps,
    mode = "both",
    increment = final_increment,
    primer_start    = study_s,
    primer_end      = study_e,
    consensus_start = safe_consensus_start,
    consensus_end   = safe_consensus_end,
    use_alignment   = USE_ALIGNMENT,
    aligned_sequences = aligned_sim_seqs
)

# --- Calculate ecological metrics at each trim step ---
message("Step 3: Calculating metrics...")

# Build OTU tables for each trim level across all taxonomic ranks
otu_tables_per_level <- analyze_all_taxonomic_levels(otu_table, trim_results$uc_files)

# Bray-Curtis dissimilarity: measures beta-diversity change at each trim step
# relative to the untrimmed baseline. This is the primary metric for the null model.
dissim_data <- calculate_dissimilarity(
  otu_tables_per_level,
  final_num_steps,
  increment = final_increment
)

# Taxonomic retention: fraction of original taxa retained at each trim step
# (complementary alpha-diversity perspective)
retention_data <- calculate_taxonomic_retention(
  otu_tables_per_level,
  final_num_steps,
  increment = final_increment
)

# Changepoint detection: identify the trim step where degradation accelerates
# (used to define the "safe trimming threshold")
threshold_data <- calculate_changepoint_thresholds(dissim_data)

# Per-taxon impact analysis: which taxa are most sensitive to trimming
impacted_taxa_data <- analyze_taxon_impact(
  otu_tables_per_level,
  num_steps = final_num_steps,
  increment = final_increment
)

#===============================================================================
# SECTION 7: OUTPUT
#===============================================================================
# Package all results into a single RDS file. The downstream aggregation script
# (07_aggregate.R) will collect results across all seeds to build the null
# distribution envelope (e.g., median + 95% CI of dissimilarity per trim step).

final_output <- list(
  task_id = TASK_ID,
  seed = SEED,
  is_outlier = is_outlier,
  num_steps = final_num_steps,
  increment = final_increment,
  step_mode = current_step_mode,
  dissim_data = dissim_data,
  retention_data = retention_data,
  threshold_data = threshold_data,
  taxon_impacts = impacted_taxa_data,
  simulation_method = "custom_grinder_inspired",
  simulation_params = list(
    n_taxa = n_taxa,
    abundance_model = if(exists("SIMULATION_ABUNDANCE_MODEL")) SIMULATION_ABUNDANCE_MODEL else "lognormal",
    pcr_bias = if(exists("SIMULATION_ADD_PCR_BIAS")) SIMULATION_ADD_PCR_BIAS else TRUE,
    pcr_cycles = if(exists("SIMULATION_PCR_CYCLES")) SIMULATION_PCR_CYCLES else 25,
    error_rate = if(exists("SIMULATION_ERROR_RATE")) SIMULATION_ERROR_RATE else 0.003,
    chimera_rate = if(exists("SIMULATION_CHIMERA_RATE")) SIMULATION_CHIMERA_RATE else 0.02,
    community_metrics = sim_community$metrics
  )
)

saveRDS(final_output, file.path(output_dir, "results.rds"))
message("Done.")
