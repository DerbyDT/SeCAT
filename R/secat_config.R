#!/usr/bin/env Rscript
# =============================================================================
# secat_config.R (Nextflow mode) — SeCAT v4.0
#
# Replaces the original secat_config.R. All parameters are read from
# environment variables set by the Nextflow process blocks (export SECAT_*).
# Variable names are identical to the original secat_config.R — no other
# R scripts need modification.
# =============================================================================

.secat_env <- function(var, default, type = "character") {
  val <- Sys.getenv(var, unset = NA_character_)
  if (is.na(val) || val == "") {
    val <- default
  } else {
    if (type == "numeric")  val <- as.numeric(val)
    if (type == "integer")  val <- as.integer(val)
    if (type == "logical")  val <- as.logical(toupper(val))
  }
  val
}

message("Loading SeCAT configuration (Nextflow mode)...")

REFERENCE_DB_PATH             <- .secat_env("SECAT_REFERENCE_DB", "")
SIMULATION_DB_PATH            <- "output/intermediate/simulation_reference_subset.fasta"
SIMULATION_USE_PREBUILT_SUBSET <- .secat_env("SECAT_SIM_USE_PREBUILT", TRUE,  "logical")
SIMULATION_MAX_SILVA_SUBSET   <- .secat_env("SECAT_SIM_MAX_SILVA_SUBSET", 10000, "integer")
VSEARCH_PATH                  <- .secat_env("SECAT_VSEARCH_PATH",  "vsearch")
SECAT_MANIFEST_PATH           <- .secat_env("SECAT_MANIFEST",      "secat_manifest_clean.tsv")
OUTDIR                        <- .secat_env("SECAT_OUTDIR",         "output")
PRIMER_DB_DIR                 <- file.path(OUTDIR, "primer_databases")
REAL_DATA_RESULTS_DIR         <- file.path(OUTDIR, "real_data_results")
SIMULATION_RESULTS_DIR        <- file.path(OUTDIR, "simulation_results")
AGGREGATED_DATA_DIR           <- file.path(OUTDIR, "aggregated_data")
FINAL_PLOTS_DIR               <- file.path(OUTDIR, "final_plots")
TAXONOMIC_LEVELS_TO_ANALYSE   <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")
NUM_SIMULATIONS_PER_PRIMER    <- .secat_env("SECAT_NUM_SIMULATIONS",       100,  "integer")
TRIM_STEP_MODE                <- .secat_env("SECAT_TRIM_STEP_MODE",         "scaled")
TRIM_INCREMENT                <- .secat_env("SECAT_TRIM_INCREMENT",          5,    "numeric")
DEFAULT_MAX_TRIM_STEPS        <- .secat_env("SECAT_DEFAULT_MAX_TRIM_STEPS",  50,   "integer")
CONSENSUS_BUFFER_STEPS        <- .secat_env("SECAT_CONSENSUS_BUFFER_STEPS",  20,   "integer")
MAX_ABSOLUTE_TRIM_STEPS       <- .secat_env("SECAT_MAX_ABSOLUTE_TRIM_STEPS", 2000, "integer")
MIN_RELATIVE_ABUNDANCE        <- .secat_env("SECAT_MIN_RELATIVE_ABUNDANCE",  0,    "numeric")
MIN_TAXA_FOR_BRAY_CURTIS      <- .secat_env("SECAT_MIN_TAXA_FOR_BRAY",       3,    "integer")
MIN_MEDIAN_READ_DEPTH         <- .secat_env("SECAT_MIN_MEDIAN_READ_DEPTH",   100,  "numeric")
TOP_N_IMPACTED_TAXA           <- 10
TOP_N_CORE_TAXA               <- 5
CORE_PREVALENCE_THRESHOLD     <- 0.90
CHANGEPOINT_PENALTY_METHOD    <- .secat_env("SECAT_CHANGEPOINT_METHOD",      "MANUAL")
CHANGEPOINT_PENALTY_MULTIPLIER <- .secat_env("SECAT_CHANGEPOINT_MULTIPLIER", 1,   "numeric")
NULL_MODEL_P_THRESHOLD        <- .secat_env("SECAT_NULL_MODEL_P",            0.05, "numeric")
NULL_MODEL_MIN_CONSECUTIVE    <- .secat_env("SECAT_NULL_MODEL_MIN_CONSEC",   3,    "integer")
NULL_MODEL_MIN_TRIM_BP        <- .secat_env("SECAT_NULL_MODEL_MIN_TRIM_BP",  5,    "numeric")
DISTANCE_CUTOFF_THRESHOLD     <- .secat_env("SECAT_DISTANCE_CUTOFF",         0.15, "numeric")
DISTANCE_CUTOFF_MIN_TRIM_BP   <- .secat_env("SECAT_DISTANCE_CUTOFF_MIN_BP",  5,    "numeric")
PLOT_DPI                      <- .secat_env("SECAT_PLOT_DPI",                300,  "integer")
FORCE_REGENERATE              <- .secat_env("SECAT_FORCE_REGENERATE",        FALSE,"logical")
SKIP_INDIVIDUAL               <- FALSE
REFERENCE_ALIGNMENT_MODE      <- .secat_env("SECAT_ALIGNMENT_MODE",          "subset")
REFERENCE_SUBSET_SIZE         <- .secat_env("SECAT_SUBSET_SIZE",             5000, "integer")
MERGE_METHOD                  <- .secat_env("SECAT_MERGE_METHOD",            "advanced")
HARMONIZE_METADATA            <- .secat_env("SECAT_HARMONIZE_METADATA",      TRUE, "logical")
METADATA_SYNONYMS <- list(
  "Latitude"   = c("lat","latitude","Lat","LAT","lat_deg"),
  "Longitude"  = c("lon","long","longitude","Lon","LONG","lon_deg"),
  "Date"       = c("date","Date","collection_date","time","Time"),
  "Depth"      = c("depth","Depth","water_depth","depth_m"),
  "Site"       = c("site","Site","location","Location","station","Station_ID"),
  "SampleType" = c("type","Type","sample_type","material","env_material")
)
CONSENSUS_OPTIMIZATION_THRESHOLD <- .secat_env("SECAT_CONSENSUS_OPT_THRESHOLD", 0.20, "numeric")
MIN_CONSENSUS_STUDIES            <- .secat_env("SECAT_MIN_CONSENSUS_STUDIES",    3,    "integer")
MIN_CONSENSUS_COVERAGE           <- .secat_env("SECAT_MIN_CONSENSUS_COVERAGE",   0.50, "numeric")
ANALYSIS_MODE                 <- .secat_env("SECAT_ANALYSIS_MODE",           "study")
STUDY_ALIGNMENT_METHOD        <- .secat_env("SECAT_ALIGNMENT_METHOD",        "modal")
USE_ALL_ASVS_FOR_MAFFT        <- .secat_env("SECAT_USE_ALL_ASVS",           TRUE,  "logical")
ASV_SAMPLE_SIZE               <- .secat_env("SECAT_ASV_SAMPLE_SIZE",        500,   "integer")
MAX_PRIMER_MISMATCH           <- .secat_env("SECAT_MAX_PRIMER_MISMATCH",    4,     "integer")
SIMULATION_ABUNDANCE_MODEL    <- .secat_env("SECAT_SIM_ABUNDANCE_MODEL",    "lognormal")
SIMULATION_LOGNORMAL_MU       <- .secat_env("SECAT_SIM_LOGNORMAL_MU",      5,     "numeric")
SIMULATION_LOGNORMAL_SIGMA    <- .secat_env("SECAT_SIM_LOGNORMAL_SIGMA",   2,     "numeric")
SIMULATION_ADD_PCR_BIAS       <- .secat_env("SECAT_SIM_ADD_PCR_BIAS",      TRUE,  "logical")
SIMULATION_PCR_CYCLES         <- .secat_env("SECAT_SIM_PCR_CYCLES",        25,    "integer")
SIMULATION_PCR_GC_BIAS_STRENGTH <- .secat_env("SECAT_SIM_PCR_GC_BIAS",    0.65,  "numeric")
SIMULATION_PCR_OPTIMAL_GC     <- .secat_env("SECAT_SIM_PCR_OPTIMAL_GC",   0.50,  "numeric")
SIMULATION_ADD_SEQUENCING_ERRORS <- .secat_env("SECAT_SIM_ADD_ERRORS",     TRUE,  "logical")
SIMULATION_ERROR_RATE         <- .secat_env("SECAT_SIM_ERROR_RATE",        0.003, "numeric")
SIMULATION_INDEL_RATE         <- .secat_env("SECAT_SIM_INDEL_RATE",        0.00003,"numeric")
SIMULATION_ERROR_POSITION_BIAS <- .secat_env("SECAT_SIM_ERROR_POSITION_BIAS", TRUE,"logical")
SIMULATION_ADD_CHIMERAS       <- .secat_env("SECAT_SIM_ADD_CHIMERAS",      FALSE, "logical")
SIMULATION_CHIMERA_RATE       <- .secat_env("SECAT_SIM_CHIMERA_RATE",      0.02,  "numeric")

message("✓ Configuration loaded successfully (Nextflow mode)")
message(strrep("=", 61))
message("  MESAP Pipeline Configuration Summary")
message(strrep("=", 61))
message(paste("  Reference Database: ", basename(REFERENCE_DB_PATH)))
message(paste("  Analysis Mode:      ", ANALYSIS_MODE))
message(paste("  VSEARCH Path:       ", VSEARCH_PATH))
message(paste("  Changepoint Method: ", CHANGEPOINT_PENALTY_METHOD))
message(paste("  Alignment Method:   ", STUDY_ALIGNMENT_METHOD))
message(paste("  Simulations/set:    ", NUM_SIMULATIONS_PER_PRIMER))
message(paste("  Trim step mode:     ", TRIM_STEP_MODE))
message(strrep("=", 61))
message("  -> Simulation mode: Custom (Grinder-inspired)")
message(paste("     PCR bias:", SIMULATION_ADD_PCR_BIAS,
              "| Errors:", SIMULATION_ADD_SEQUENCING_ERRORS,
              "| Chimeras:", SIMULATION_ADD_CHIMERAS))
