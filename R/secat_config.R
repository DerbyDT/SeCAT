#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   secat_config.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# PURPOSE:  Central configuration hub — reads all pipeline parameters from
#           Nextflow-exported environment variables into R global variables.
#
# OVERVIEW:
#   This script replaces hard-coded R configuration with a Nextflow-driven
#   approach. Every parameter is read from SECAT_* environment variables that
#   Nextflow exports into each process block. The variable names are kept
#   identical to the original secat_config.R so that downstream R scripts
#   (simulation, trimming, changepoint, consensus, merge) require zero changes.
#
# USAGE:
#   Sourced automatically by each R process in main.nf via:
#     source("R/secat_config.R")
#   Not intended to be run interactively — parameters come from params.yaml
#   via Nextflow environment variable injection.
# ==============================================================================


# ------------------------------------------------------------------------------
# Helper: .secat_env()
# Reads a single environment variable exported by Nextflow, falling back to a
# default value if unset or empty. Coerces to the requested type.
#
# Arguments:
#   var     — environment variable name (string, e.g. "SECAT_NUM_SIMULATIONS")
#   default — fallback value when the env var is absent
#   type    — one of "character", "numeric", "integer", "logical"
# ------------------------------------------------------------------------------
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

# ==============================================================================
# SECTION 1: PATHS AND DIRECTORIES
# ==============================================================================

# Path to the SILVA reference alignment (full-length, aligned FASTA).
# Valid values: absolute path to SILVA_138.x_SSURef_NR99 aligned FASTA.
# Data impact: this defines the coordinate space for all trimming decisions.
#   Using a different SILVA version changes position numbering, invalidating
#   any previously computed consensus regions.
REFERENCE_DB_PATH             <- .secat_env("SECAT_REFERENCE_DB", "")

# Path to the subset of SILVA used for simulation. Generated automatically
# during the SUBSET_SILVA process. Fixed relative path — not user-configurable.
SIMULATION_DB_PATH            <- "output/intermediate/simulation_reference_subset.fasta"

# Whether to reuse a previously built SILVA subset for simulations.
# Valid values: TRUE / FALSE.
# Data impact: TRUE skips re-subsetting (faster reruns); FALSE forces fresh
#   random sampling, which may give slightly different null distributions.
SIMULATION_USE_PREBUILT_SUBSET <- .secat_env("SECAT_SIM_USE_PREBUILT", TRUE,  "logical")

# Maximum number of SILVA sequences to retain in the simulation reference subset.
# Valid values: positive integer (typically 5000-50000).
# Data impact: larger subsets produce more realistic null communities but
#   increase memory use and runtime linearly. Current default (10000) balances
#   taxonomic coverage against compute cost.
SIMULATION_MAX_SILVA_SUBSET   <- .secat_env("SECAT_SIM_MAX_SILVA_SUBSET", 10000, "integer")

# Path to the VSEARCH binary for chimera detection and sequence searching.
# Valid values: "vsearch" (if on PATH) or an absolute path to the binary.
VSEARCH_PATH                  <- .secat_env("SECAT_VSEARCH_PATH",  "vsearch")

# Path to the SeCAT manifest TSV listing studies, primer pairs, and ASV tables.
# Valid values: absolute or relative path to a correctly formatted manifest.
# Data impact: this is the primary input — defines which studies enter the pipeline.
SECAT_MANIFEST_PATH           <- .secat_env("SECAT_MANIFEST",      "secat_manifest_clean.tsv")

# Root output directory. All subdirectories are created beneath this.
OUTDIR                        <- .secat_env("SECAT_OUTDIR",         "output")

# Derived output subdirectories (not user-configurable).
PRIMER_DB_DIR                 <- file.path(OUTDIR, "primer_databases")
REAL_DATA_RESULTS_DIR         <- file.path(OUTDIR, "real_data_results")
SIMULATION_RESULTS_DIR        <- file.path(OUTDIR, "simulation_results")
AGGREGATED_DATA_DIR           <- file.path(OUTDIR, "aggregated_data")
FINAL_PLOTS_DIR               <- file.path(OUTDIR, "final_plots")

# ==============================================================================
# SECTION 2: TAXONOMIC ANALYSIS SETTINGS
# ==============================================================================

# Taxonomic ranks at which trimming impact is assessed via Bray-Curtis distance.
# Each rank produces independent degradation curves and changepoint tests.
# Data impact: including "ASV" captures fine-scale compositional shifts but is
#   noisier; coarser ranks (Phylum, Class) are more stable but less sensitive.
TAXONOMIC_LEVELS_TO_ANALYSE   <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

# ==============================================================================
# SECTION 3: SIMULATION ENGINE PARAMETERS
# Generates synthetic communities with Grinder-inspired PCR/sequencing artefacts
# to build null distributions for each primer pair's trimming response.
# ==============================================================================

# Number of independent synthetic communities per primer pair.
# Valid values: positive integer (50-1000; 100 recommended for publication).
# Data impact: more simulations tighten the null distribution and reduce
#   stochastic variability in p-values. Fewer than 50 may give unreliable
#   Type I error control. Runtime scales linearly. Default calibrated for
#   optimal Type I error control via testing.
NUM_SIMULATIONS_PER_PRIMER    <- .secat_env("SECAT_NUM_SIMULATIONS",       100,  "integer")

# ==============================================================================
# SECTION 4: TRIMMING CONFIGURATION
# Controls how progressively more aggressive trimming is applied to aligned ASVs.
# ==============================================================================

# Trim step mode: how trim positions are spaced.
# Valid values: "scaled" (steps proportional to amplicon length) or
#   "absolute" (fixed base-pair increments).
# Data impact: "scaled" ensures comparable trimming across primers of different
#   lengths; "absolute" gives uniform bp resolution but may over/under-sample
#   short/long amplicons. Default calibrated via testing.
TRIM_STEP_MODE                <- .secat_env("SECAT_TRIM_STEP_MODE",         "scaled")

# Base-pair increment per trim step (used in "absolute" mode, or as the base
# increment in "scaled" mode).
# Valid values: positive numeric (typically 5-100 bp).
# Data impact: smaller increments give finer resolution for detecting the
#   degradation threshold but increase the number of steps (and compute time).
TRIM_INCREMENT                <- .secat_env("SECAT_TRIM_INCREMENT",          5,    "numeric")

# Maximum number of trim steps to evaluate per primer/study.
# Valid values: positive integer (10-200).
# Data impact: caps the trimming search space. Too few may miss the degradation
#   point; too many wastes compute on irrelevant extreme trims.
DEFAULT_MAX_TRIM_STEPS        <- .secat_env("SECAT_DEFAULT_MAX_TRIM_STEPS",  50,   "integer")

# Extra trim steps beyond the consensus region boundary for safety margin.
# Valid values: positive integer (5-50).
# Data impact: ensures degradation curves extend past the candidate trim point,
#   giving changepoint detection sufficient post-change data.
CONSENSUS_BUFFER_STEPS        <- .secat_env("SECAT_CONSENSUS_BUFFER_STEPS",  20,   "integer")

# Hard ceiling on total trim steps (guards against runaway in edge cases).
# Valid values: positive integer.
# Data impact: pure safety limit; should never be reached in normal operation.
MAX_ABSOLUTE_TRIM_STEPS       <- .secat_env("SECAT_MAX_ABSOLUTE_TRIM_STEPS", 2000, "integer")

# ==============================================================================
# SECTION 5: QUALITY FILTERS
# Pre-trim filters that remove low-quality or uninformative samples/taxa.
# ==============================================================================

# Minimum relative abundance for a taxon to be retained before trimming analysis.
# Valid values: 0-1 (0 = no filtering; 0.001 = drop taxa below 0.1%).
# Data impact: higher values remove rare taxa, reducing noise but potentially
#   masking impacts on low-abundance community members. 0 retains everything.
MIN_RELATIVE_ABUNDANCE        <- .secat_env("SECAT_MIN_RELATIVE_ABUNDANCE",  0,    "numeric")

# Minimum number of taxa required to compute Bray-Curtis distance.
# Valid values: positive integer (>= 2, typically 3-5).
# Data impact: samples with fewer taxa than this are excluded from degradation
#   curves. Prevents degenerate distance calculations from near-empty communities.
MIN_TAXA_FOR_BRAY_CURTIS      <- .secat_env("SECAT_MIN_TAXA_FOR_BRAY",       3,    "integer")

# Minimum median read depth across samples for a study to be included.
# Valid values: positive numeric (typically 100-1000).
# Data impact: filters out severely undersequenced studies that would produce
#   unreliable compositional profiles. Lower thresholds include more studies
#   but risk noisy distance estimates.
MIN_MEDIAN_READ_DEPTH         <- .secat_env("SECAT_MIN_MEDIAN_READ_DEPTH",   100,  "numeric")

# ==============================================================================
# SECTION 6: REPORTING PARAMETERS (hard-coded, not Nextflow-configurable)
# ==============================================================================

# Number of most-impacted taxa to highlight in per-study impact reports.
TOP_N_IMPACTED_TAXA           <- 10

# Number of core taxa to report per study.
TOP_N_CORE_TAXA               <- 5

# Prevalence threshold defining "core" taxa (fraction of samples where present).
# Valid values: 0-1 (0.90 = taxon must appear in >= 90% of samples).
CORE_PREVALENCE_THRESHOLD     <- 0.90

# ==============================================================================
# SECTION 7: CHANGEPOINT DETECTION (PELT ALGORITHM)
# Identifies the trim step where Bray-Curtis distance shows a statistically
# significant shift — the point beyond which trimming causes unacceptable
# community distortion.
# ==============================================================================

# Penalty method for the PELT (Pruned Exact Linear Time) changepoint algorithm.
# Valid values: "MANUAL" (user-scaled penalty), "BIC", "AIC", "MBIC", "SIC".
# Data impact: "MANUAL" with the multiplier gives direct control over
#   sensitivity. Information-criterion methods (BIC, AIC) auto-calibrate but
#   may be too liberal or conservative depending on the noise structure.
#   Default calibrated via testing for optimal Type I error control.
CHANGEPOINT_PENALTY_METHOD    <- .secat_env("SECAT_CHANGEPOINT_METHOD",      "MANUAL")

# Multiplier applied to the MANUAL penalty (penalty = multiplier * var(data)).
# Valid values: positive numeric (0.1-100; 1 = baseline).
# Data impact: increasing makes PELT less sensitive (fewer changepoints,
#   tolerates more degradation before flagging); decreasing makes it more
#   sensitive (flags smaller shifts, higher false positive rate).
#   Default calibrated via testing for optimal Type I error control.
CHANGEPOINT_PENALTY_MULTIPLIER <- .secat_env("SECAT_CHANGEPOINT_MULTIPLIER", 1,   "numeric")

# ==============================================================================
# SECTION 8: NULL MODEL COMPARISON
# Tests whether observed trimming degradation exceeds what synthetic (null)
# communities show — i.e., whether the real signal is biologically meaningful
# or just stochastic noise.
# ==============================================================================

# P-value threshold for rejecting the null hypothesis that observed degradation
# equals simulated degradation.
# Valid values: 0-1 (0.05 = standard significance level).
# Data impact: lower values are more conservative (fewer studies flagged as
#   significantly degraded), reducing false positives but potentially missing
#   real trimming artefacts.
#   Default calibrated via testing for optimal Type I error control.
NULL_MODEL_P_THRESHOLD        <- .secat_env("SECAT_NULL_MODEL_P",            0.05, "numeric")

# Number of consecutive trim steps that must exceed the null threshold to
# trigger a "significant degradation" verdict.
# Valid values: positive integer (1-10; 3 recommended).
# Data impact: higher values require sustained degradation (not a single noisy
#   spike), reducing false positives but delaying detection.
NULL_MODEL_MIN_CONSECUTIVE    <- .secat_env("SECAT_NULL_MODEL_MIN_CONSEC",   3,    "integer")

# Minimum base pairs of trimming before null model testing begins.
# Valid values: positive numeric (typically 5-20 bp).
# Data impact: prevents flagging trivially small trims where distance changes
#   are unreliable. Acts as a "dead zone" at the start of the trimming curve.
NULL_MODEL_MIN_TRIM_BP        <- .secat_env("SECAT_NULL_MODEL_MIN_TRIM_BP",  5,    "numeric")

# ==============================================================================
# SECTION 9: DISTANCE CUTOFF METHOD (SECONDARY THRESHOLD)
# An absolute Bray-Curtis distance threshold — a simpler complement to the
# null-model approach. A study is flagged if median Bray-Curtis exceeds this
# value at any trim step.
# ==============================================================================

# Absolute Bray-Curtis distance threshold for the cutoff method.
# Valid values: 0-1 (0.15 = 15% compositional dissimilarity).
# Data impact: lower values are stricter (flag more studies); higher values
#   tolerate more compositional change before trimming is deemed harmful.
#   0.15 represents a moderate ecological threshold where community structure
#   shifts are likely biologically meaningful.
DISTANCE_CUTOFF_THRESHOLD     <- .secat_env("SECAT_DISTANCE_CUTOFF",         0.15, "numeric")

# Minimum base pairs before the distance cutoff is evaluated.
# Valid values: positive numeric.
# Data impact: same dead-zone logic as NULL_MODEL_MIN_TRIM_BP.
DISTANCE_CUTOFF_MIN_TRIM_BP   <- .secat_env("SECAT_DISTANCE_CUTOFF_MIN_BP",  5,    "numeric")

# ==============================================================================
# SECTION 10: OUTPUT AND REGENERATION CONTROL
# ==============================================================================

# Plot resolution in dots per inch. Affects file size and print quality.
# Valid values: positive integer (150-600; 300 = publication standard).
PLOT_DPI                      <- .secat_env("SECAT_PLOT_DPI",                300,  "integer")

# Force regeneration of all outputs, even if cached results exist.
# Valid values: TRUE / FALSE.
# Data impact: TRUE ensures a fully fresh run (useful after parameter changes);
#   FALSE reuses existing outputs where possible (faster for iterative runs).
FORCE_REGENERATE              <- .secat_env("SECAT_FORCE_REGENERATE",        FALSE,"logical")

# Skip individual study-level analysis (used internally for aggregation-only runs).
SKIP_INDIVIDUAL               <- FALSE

# ==============================================================================
# SECTION 11: REFERENCE ALIGNMENT SETTINGS
# Controls how ASVs are aligned to the SILVA reference for coordinate mapping.
# ==============================================================================

# Alignment mode: how much of the SILVA reference to use for DECIPHER alignment.
# Valid values: "subset" (random SILVA subset for speed) or "full" (entire DB).
# Data impact: "subset" is dramatically faster and usually sufficient; "full"
#   may improve alignment accuracy for divergent sequences but requires
#   substantially more memory (64 GB+) and time.
REFERENCE_ALIGNMENT_MODE      <- .secat_env("SECAT_ALIGNMENT_MODE",          "subset")

# Number of SILVA sequences in the alignment reference subset.
# Valid values: positive integer (1000-50000; 5000 recommended).
# Data impact: larger subsets improve alignment accuracy at the cost of memory
#   and time. Diminishing returns above ~10000 for typical marine/soil studies.
REFERENCE_SUBSET_SIZE         <- .secat_env("SECAT_SUBSET_SIZE",             5000, "integer")

# ==============================================================================
# SECTION 12: MERGE AND METADATA HARMONISATION
# Controls how trimmed studies are combined into the final unified dataset.
# ==============================================================================

# Merge method for combining trimmed ASV tables across studies.
# Valid values: "advanced" (taxonomy-aware merging with conflict resolution)
#   or "simple" (direct concatenation).
# Data impact: "advanced" reconciles taxonomic assignments across studies and
#   handles ASV identity conflicts; "simple" is faster but may introduce
#   duplicated taxa with inconsistent annotations.
MERGE_METHOD                  <- .secat_env("SECAT_MERGE_METHOD",            "advanced")

# Whether to standardise metadata column names across studies.
# Valid values: TRUE / FALSE.
# Data impact: TRUE applies the synonym mapping below to unify column names
#   (e.g., "lat" -> "Latitude"), enabling cross-study metadata joins.
HARMONIZE_METADATA            <- .secat_env("SECAT_HARMONIZE_METADATA",      TRUE, "logical")

# Metadata column name synonyms. Each canonical name (key) maps to a list of
# common variants found in published 16S metadata. Used during harmonisation
# to standardise column names across heterogeneous study metadata.
METADATA_SYNONYMS <- list(
  "Latitude"   = c("lat","latitude","Lat","LAT","lat_deg"),
  "Longitude"  = c("lon","long","longitude","Lon","LONG","lon_deg"),
  "Date"       = c("date","Date","collection_date","time","Time"),
  "Depth"      = c("depth","Depth","water_depth","depth_m"),
  "Site"       = c("site","Site","location","Location","station","Station_ID"),
  "SampleType" = c("type","Type","sample_type","material","env_material")
)

# ==============================================================================
# SECTION 13: CONSENSUS REGION OPTIMISATION
# Finds the maximal overlapping 16S region across all primer pairs, subject to
# a maximum tolerable degradation threshold.
# ==============================================================================

# Maximum tolerable Bray-Curtis degradation for including a study in the
# consensus trim region.
# Valid values: 0-1 (0.20 = allow up to 20% compositional change).
# Data impact: higher values include more studies in the consensus (wider
#   participation) but accept greater distortion. Lower values produce a more
#   conservative consensus that preserves ecological signal but may exclude
#   studies with narrow amplicons.
CONSENSUS_OPTIMIZATION_THRESHOLD <- .secat_env("SECAT_CONSENSUS_OPT_THRESHOLD", 0.20, "numeric")

# Minimum number of studies that must contribute to the consensus region.
# Valid values: positive integer (>= 2).
# Data impact: higher values ensure the consensus is broadly supported but
#   may fail if too few studies share overlapping regions.
MIN_CONSENSUS_STUDIES            <- .secat_env("SECAT_MIN_CONSENSUS_STUDIES",    3,    "integer")

# Minimum fraction of studies that must be representable in the consensus.
# Valid values: 0-1 (0.50 = at least half of input studies).
# Data impact: acts as a proportional guard — ensures the consensus region is
#   not driven by a small minority of studies.
MIN_CONSENSUS_COVERAGE           <- .secat_env("SECAT_MIN_CONSENSUS_COVERAGE",   0.50, "numeric")

# ==============================================================================
# SECTION 14: STUDY-LEVEL ALIGNMENT AND ASV HANDLING
# ==============================================================================

# Analysis granularity: whether trimming decisions are made per-study or
# per-primer-pair.
# Valid values: "study" (each study assessed independently) or "primer"
#   (all studies sharing a primer pair assessed together).
# Data impact: "study" mode captures study-specific biases (e.g., different
#   extraction protocols); "primer" mode pools data for more statistical power
#   but may mask study-level effects.
ANALYSIS_MODE                 <- .secat_env("SECAT_ANALYSIS_MODE",           "study")

# Method for determining the representative aligned position per ASV per study.
# Valid values: "modal" (most common alignment position) or "median".
# Data impact: "modal" is robust to outlier alignments; "median" gives a
#   central tendency that may be pulled by skewed distributions.
STUDY_ALIGNMENT_METHOD        <- .secat_env("SECAT_ALIGNMENT_METHOD",        "modal")

# Whether to use all ASVs for DECIPHER alignment or subsample.
# Valid values: TRUE / FALSE.
# Data impact: TRUE uses every ASV (most accurate but memory-intensive);
#   FALSE subsamples to ASV_SAMPLE_SIZE (faster, lower memory, slight
#   alignment accuracy trade-off for very diverse datasets).
USE_ALL_ASVS                  <- .secat_env("SECAT_USE_ALL_ASVS",           TRUE,  "logical")

# Number of ASVs to subsample when USE_ALL_ASVS is FALSE.
# Valid values: positive integer (100-10000).
# Data impact: more ASVs give better alignment coverage; fewer are faster.
#   Only relevant when USE_ALL_ASVS = FALSE.
ASV_SAMPLE_SIZE               <- .secat_env("SECAT_ASV_SAMPLE_SIZE",        500,   "integer")

# Maximum allowed mismatches when matching primer sequences to study ASVs.
# Valid values: non-negative integer (0-6; 4 recommended).
# Data impact: higher values tolerate more primer degeneracy / sequencing
#   errors in the primer-binding region, including more ASVs but risking
#   false primer assignments.
MAX_PRIMER_MISMATCH           <- .secat_env("SECAT_MAX_PRIMER_MISMATCH",    4,     "integer")

# ==============================================================================
# SECTION 15: SIMULATION ENGINE — COMMUNITY AND ERROR MODEL
# These parameters control the Grinder-inspired synthetic community generator.
# Together they define the null expectation: what trimming degradation looks
# like when there is no real biological signal to lose.
# ==============================================================================

# Species abundance distribution model for synthetic communities.
# Valid values: "lognormal" (realistic SAD) or "uniform" (equal abundances).
# Data impact: "lognormal" produces communities dominated by a few abundant
#   taxa (ecologically realistic); "uniform" gives equal weight to all taxa,
#   which may underestimate trimming impacts on rare-taxon-driven distance.
SIMULATION_ABUNDANCE_MODEL    <- .secat_env("SECAT_SIM_ABUNDANCE_MODEL",    "lognormal")

# Mean of the log-normal species abundance distribution (log scale).
# Valid values: positive numeric (typically 3-8).
# Data impact: higher values shift the SAD toward higher mean abundances,
#   producing more even communities. Lower values create stronger dominance
#   by a few taxa. Default (5) mimics typical marine/soil 16S profiles.
SIMULATION_LOGNORMAL_MU       <- .secat_env("SECAT_SIM_LOGNORMAL_MU",      5,     "numeric")

# Standard deviation of the log-normal SAD (log scale).
# Valid values: positive numeric (typically 1-4).
# Data impact: higher values increase abundance evenness spread, creating
#   communities with a wider range of rare-to-dominant taxa. Default (2)
#   produces realistic dominance hierarchies.
SIMULATION_LOGNORMAL_SIGMA    <- .secat_env("SECAT_SIM_LOGNORMAL_SIGMA",   2,     "numeric")

# Whether to simulate GC-content-dependent PCR amplification bias.
# Valid values: TRUE / FALSE.
# Data impact: TRUE introduces realistic compositional distortion from
#   differential amplification efficiency. FALSE assumes perfect PCR
#   (unrealistic but useful for isolating other effects).
SIMULATION_ADD_PCR_BIAS       <- .secat_env("SECAT_SIM_ADD_PCR_BIAS",      TRUE,  "logical")

# Number of PCR cycles to simulate.
# Valid values: positive integer (15-35; 25 = standard for most 16S protocols).
# Data impact: more cycles amplify GC bias exponentially. Values above 30
#   produce extreme compositional skew.
SIMULATION_PCR_CYCLES         <- .secat_env("SECAT_SIM_PCR_CYCLES",        25,    "integer")

# Strength of GC-content bias during simulated PCR.
# Valid values: 0-1 (0 = no bias; 1 = extreme bias; 0.65 = moderate).
# Data impact: higher values cause stronger differential amplification of
#   sequences near the optimal GC content, distorting community composition.
#   Default (0.65) reflects empirically measured PCR bias in Illumina 16S
#   library preparation.
SIMULATION_PCR_GC_BIAS_STRENGTH <- .secat_env("SECAT_SIM_PCR_GC_BIAS",    0.65,  "numeric")

# Optimal GC content for simulated PCR efficiency (fraction).
# Valid values: 0-1 (0.50 = balanced; typical for most polymerases).
# Data impact: sequences with GC content near this value are preferentially
#   amplified. Shifting away from 0.50 biases toward AT-rich or GC-rich taxa.
SIMULATION_PCR_OPTIMAL_GC     <- .secat_env("SECAT_SIM_PCR_OPTIMAL_GC",   0.50,  "numeric")

# Whether to add simulated sequencing errors (substitutions and indels).
# Valid values: TRUE / FALSE.
# Data impact: TRUE produces realistic noisy reads (may inflate ASV counts
#   and blur distances); FALSE gives error-free sequences.
SIMULATION_ADD_SEQUENCING_ERRORS <- .secat_env("SECAT_SIM_ADD_ERRORS",     TRUE,  "logical")

# Per-base substitution error rate for simulated sequencing.
# Valid values: 0-1 (0.003 = 0.3%, typical Illumina MiSeq error rate).
# Data impact: higher rates increase noise in simulated communities, widening
#   null distributions and making it harder to detect real degradation signals.
SIMULATION_ERROR_RATE         <- .secat_env("SECAT_SIM_ERROR_RATE",        0.003, "numeric")

# Per-base insertion/deletion error rate.
# Valid values: 0-1 (0.00003 = 0.003%, typical for Illumina platforms).
# Data impact: indels cause frameshift-like artefacts in alignments; higher
#   rates produce more spurious ASVs in the null model.
SIMULATION_INDEL_RATE         <- .secat_env("SECAT_SIM_INDEL_RATE",        0.00003,"numeric")

# Whether sequencing error rate increases toward the 3' end of reads.
# Valid values: TRUE / FALSE.
# Data impact: TRUE mimics real Illumina quality decay, concentrating errors
#   at read ends (which are the positions most affected by trimming).
#   Produces more realistic null distributions for trimming impact assessment.
SIMULATION_ERROR_POSITION_BIAS <- .secat_env("SECAT_SIM_ERROR_POSITION_BIAS", TRUE,"logical")

# Whether to add simulated chimeric sequences.
# Valid values: TRUE / FALSE.
# Data impact: TRUE introduces chimeras (inter-template PCR artefacts),
#   inflating diversity in the null model. FALSE assumes perfect chimera
#   removal by upstream tools (e.g., DADA2). Disabled by default because
#   most modern pipelines remove chimeras before SeCAT input.
SIMULATION_ADD_CHIMERAS       <- .secat_env("SECAT_SIM_ADD_CHIMERAS",      FALSE, "logical")

# Fraction of reads that are chimeric when chimera simulation is enabled.
# Valid values: 0-1 (0.02 = 2%, conservative estimate).
# Data impact: only active when SIMULATION_ADD_CHIMERAS = TRUE. Higher rates
#   increase spurious diversity in null communities.
SIMULATION_CHIMERA_RATE       <- .secat_env("SECAT_SIM_CHIMERA_RATE",      0.02,  "numeric")

# ==============================================================================
# CONFIGURATION SUMMARY (printed to Nextflow process log)
# ==============================================================================
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
