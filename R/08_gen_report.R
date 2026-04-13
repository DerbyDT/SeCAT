#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   scripts/08_gen_report.R
# PIPELINE: SeCAT (Sequence Consensus Analysis Tool)
# PHASE:    Phase 6: Final Report Generation
# VERSION:  2.0 (Dynamic Step Calculation)
# AUTHOR:   [Author Name]
#
# PURPOSE:
#   This is the visualization engine of SeCAT. It takes the "Verdicts" from
#   Phase 5 and the raw data from Phases 3 & 4 to generate:
#   1. Individual PDF Reports for every study (Dissimilarity, Retention, Impact).
#   2. A Master Summary PDF comparing all studies (Alignment Map, Forest Plot).
#
# CHUNK 1: SETUP & DATA LOADING
#   - Loads necessary plotting libraries (ggplot2, patchwork, magick).
#   - Loads the aggregated simulation baselines (The "Gray Ribbon").
#   - Loads the master verdict table (The "Decisions").
#   - Defines shared helper functions for plotting thresholds.
# ==============================================================================

message("--- Loading libraries and configuration ---")
suppressPackageStartupMessages({
    suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
    library(here)
    library(patchwork)   # For combining plots (plot_layout)
    library(ggrepel)     # For non-overlapping text labels
    library(scales)      # For pretty axis formatting
    library(gridExtra)   # Legacy grid tools
    library(vegan)       # Biological diversity stats
    library(grid)
    library(gtable)
    library(magick)      # For stitching PNGs into PDFs
})

# --- CRITICAL: Define log_and_flush function ---
# Ensures logs appear in real-time in the SGE output file
log_and_flush <- function(message) {
  cat(paste(Sys.time(), "|", message, "\n"))
  flush.console()
}

# ==============================================================================
# SECTION 1: CONFIGURATION
# ==============================================================================

# Load from config with defaults (allows running interactively or via script)
FORCE_REGENERATE <- if (exists("FORCE_REGENERATE")) FORCE_REGENERATE else FALSE
SKIP_INDIVIDUAL <- if (exists("SKIP_INDIVIDUAL")) SKIP_INDIVIDUAL else FALSE
PLOT_WIDTH <- if (exists("PLOT_WIDTH")) PLOT_WIDTH else 16
PLOT_HEIGHT <- if (exists("PLOT_HEIGHT")) PLOT_HEIGHT else 11
PLOT_DPI <- if (exists("PLOT_DPI")) PLOT_DPI else 300

# Load project configuration
if (file.exists(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R"))) {
    source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R"))
} else {
    # Fallback paths if not running from root
    AGGREGATED_DATA_DIR <- file.path(OUTDIR, "aggregated_data")
    SECAT_MANIFEST_PATH <- SECAT_MANIFEST_PATH
    FINAL_PLOTS_DIR <- file.path(OUTDIR, "final_plots")
}

output_dir <- if (exists("FINAL_PLOTS_DIR")) FINAL_PLOTS_DIR else file.path(OUTDIR, "final_plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# SECTION 2: LOAD SHARED DATA
# ==============================================================================
log_and_flush("--- Loading shared data ---")

# --- 2a. Load Simulation Baseline Statistics ---
# This data creates the "Blue/Gray Ribbon" representing the Null Model.
sim_baselines_path <- here(AGGREGATED_DATA_DIR, "simulation_baseline_statistics.csv")
sim_baselines <- if (file.exists(sim_baselines_path)) {
    message("✅ Found simulation baseline statistics")
    read_csv(sim_baselines_path, show_col_types = FALSE)
} else {
    message("⚠️ Simulation baseline statistics file not found - creating empty structure")
    tibble::tibble(
        task_id = character(),
        Level = character(),
        Trim_BP = numeric(),
        MeanDissimilarity = numeric(),
        CILower = numeric(),
        CIUpper = numeric()
    )
}

# --- 2b. Load Simulation Retention Curves ---
# This data allows us to plot the expected loss of taxa (Alpha Diversity).
ret_curves_path <- here(AGGREGATED_DATA_DIR, "simulation_retention_curves.csv")
sim_retention_curves_raw <- if (file.exists(ret_curves_path)) {
  message("✅ Found simulation retention curves (raw data)")
  readr::read_csv(ret_curves_path, show_col_types = FALSE)
} else {
  message("⚠️  Simulation retention curves file not found - creating empty structure")
  tibble::tibble(
    task_id = character(),
    Level = character(),
    Trim_BP = numeric(),
    Retention = numeric(),
    mode = character()
  )
}

# --- 2c. Process Retention Curves ---
# Logic: Checks if the loaded data is Raw (per-seed) or Pre-Aggregated.
# If Raw, it calculates Mean/SD/CI on the fly.
sim_retention_curves <- if (nrow(sim_retention_curves_raw) > 0) {
  message("   → Processing retention data for baseline plotting...")

  df <- sim_retention_curves_raw

  # Standardize column names (Task ID mapping)
  if ("study_name" %in% names(df) && !"task_id" %in% names(df)) {
    df <- df %>% rename(task_id = study_name)
  } else if ("primer_name" %in% names(df) && !"task_id" %in% names(df)) {
    df <- df %>% rename(task_id = primer_name)
  }

  if ("MeanRetention" %in% names(df)) {
      message("   ℹ️ Data appears to be pre-aggregated (MeanRetention column found). Using as-is.")
      # Ensure optional columns exist for compatibility
      if (!"SDRetention" %in% names(df)) df$SDRetention <- NA
      if (!"n_simulations" %in% names(df)) df$n_simulations <- NA

      result <- df %>%
          select(any_of(c("task_id", "Level", "Trim_BP", "MeanRetention", "SDRetention", "CILower", "CIUpper", "n_simulations")))

  } else {
      # Fallback: Aggregating Raw Data
      message("   ℹ️ Data appears to be raw. Aggregating now...")
      result <- df %>%
        { if ("mode" %in% names(.)) filter(., mode == "both") else . } %>%
        group_by(task_id, Level, Trim_BP) %>%
        summarise(
          n_simulations = n(),
          MeanRetention = mean(Retention, na.rm = TRUE),
          SDRetention = sd(Retention, na.rm = TRUE),
          CILower = quantile(Retention, 0.025, na.rm = TRUE),
          CIUpper = quantile(Retention, 0.975, na.rm = TRUE),
          .groups = "drop"
        )
      message(paste("   ✅ Aggregated", nrow(df), "raw rows into", nrow(result), "summary rows"))
  }
  result
} else {
  tibble::tibble()
}

# --- 2d. Load Master Verdict Table ---
# This contains the statistical results (Pass/Fail, Thresholds) for all studies.
master_verdicts_path <- here(AGGREGATED_DATA_DIR, "master_verdict_table.csv")
master_verdicts <- if (file.exists(master_verdicts_path)) {
    message("✅ Found master verdict table")
    read_csv(master_verdicts_path, show_col_types = FALSE)
} else {
    message("⚠️ Master verdict table not found")
    tibble()
}

if (nrow(master_verdicts) == 0) {
  stop("⛔ FATAL: Master verdict table is empty. Cannot generate reports.")
}

message(paste("   → Loaded", nrow(master_verdicts), "verdict records for",
              length(unique(master_verdicts$Study)), "studies"))

# --- 2e. Load Manifest ---
master_manifest <- read_tsv(SECAT_MANIFEST_PATH, show_col_types = FALSE)
master_manifest$study_name <- trimws(master_manifest$study_name)

# Clean up memory after large loads
gc()

# ==============================================================================
# SECTION 3: SHARED HELPER FUNCTIONS
# ==============================================================================

# Function: is_valid_df
# Description: Checks if an object is a non-empty data frame.
is_valid_df <- function(df) {
  !is.null(df) && inherits(df, "data.frame") && nrow(df) > 0
}

# Function: compute_first_threshold
# Description: Extracts the correct degradation threshold for summary functions,
#              using Confirmed_Threshold_BP when available (v3.2+) and falling back
#              to min() for legacy verdict tables. Used by all four summary functions
#              to replace the old inline pmap_dbl(min()) pattern.
#
#              For WARNING_SINGLE: returns NA (lone trigger, not confirmed).
#              For CONFIRMED:      returns Confirmed_Threshold_BP.
#              For NONE / legacy:  returns NA or min() respectively.
compute_first_threshold <- function(cp, cut, null, confirmed_bp = NA_real_, consensus_status = NA_character_) {
  # v3.2+ path
  if (!is.na(consensus_status)) {
    if (consensus_status %in% c("CONFIRMED", "HARD_FAIL")) return(confirmed_bp)
    return(NA_real_)   # WARNING_SINGLE or NONE
  }
  # Legacy fallback
  thresholds <- c(cp, cut, null)
  thresholds <- thresholds[!is.na(thresholds)]
  if (length(thresholds) == 0) NA_real_ else min(thresholds)
}

# Function: add_threshold_lines
# Description: Adds the visual "Goal Post" (Green) and "Failure Point" (Red) lines.
#              Critically, it divides by 'increment' to plot correctly on the X-axis,
#              which is usually in "Steps", not "BP".
add_threshold_lines <- function(p, required_thresh_bp, observed_thresh_bp,
                                increment = NULL, warning_thresh_bp = NA_real_) {
  if (is.null(increment) || is.na(increment)) {
    warning("No increment provided to add_threshold_lines(). Using default of 10 bp. Threshold lines may be misplaced.")
    increment <- 10
  }

  # Green Dashed Line: The Minimum Required Trimming (Target)
  if (!is.na(required_thresh_bp)) {
    p <- p + geom_vline(xintercept = required_thresh_bp / increment,
                        linetype = "dashed", color = "darkgreen", linewidth = 1)
  }
  # Red Dotted Line: Confirmed multi-method degradation threshold
  if (!is.na(observed_thresh_bp)) {
    p <- p + geom_vline(xintercept = observed_thresh_bp / increment,
                        linetype = "dotted", color = "darkred", linewidth = 1)
  }
  # Orange Dotted Line: Single-method WARNING_SINGLE trigger (study still passes)
  if (!is.na(warning_thresh_bp)) {
    p <- p + geom_vline(xintercept = warning_thresh_bp / increment,
                        linetype = "dotted", color = "darkorange", linewidth = 1.2)
  }
  return(p)
}

# Function: get_first_degradation_threshold
# Description: Returns the CONFIRMED degradation threshold from consensus voting
#              (v3.2+: uses Confirmed_Threshold_BP column from 07_aggregate.R).
#
#              Falls back to old min() behaviour only if the new columns are
#              absent — this keeps the script compatible with older verdict tables.
#
#              For WARNING_SINGLE status (single method triggered, no consensus):
#              returns NA so that NO red line is drawn. The warning is shown in
#              the subtitle instead. This prevents false-positive red lines.
get_first_degradation_threshold <- function(vdf) {
  if (is.null(vdf) || nrow(vdf) == 0) return(NA_real_)

  # v3.2+ path: use pre-computed consensus result
  if ("Confirmed_Threshold_BP" %in% names(vdf) && "Consensus_Status" %in% names(vdf)) {
    status <- vdf$Consensus_Status[1]
    bp     <- vdf$Confirmed_Threshold_BP[1]

    if (!is.na(status) && status %in% c("CONFIRMED", "HARD_FAIL")) return(bp)
    # WARNING_SINGLE or NONE: no confirmed threshold — return NA (no red line)
    return(NA_real_)
  }

  # Legacy fallback: old verdict tables without consensus columns
  t_cp   <- if ("Threshold_Observed_Changepoint" %in% names(vdf)) vdf$Threshold_Observed_Changepoint[1] else NA_real_
  t_cut  <- if ("Threshold_Observed_Cutoff"      %in% names(vdf)) vdf$Threshold_Observed_Cutoff[1]      else NA_real_
  t_null <- if ("Threshold_Observed_NullModel"   %in% names(vdf)) vdf$Threshold_Observed_NullModel[1]   else NA_real_
  tt <- c(t_cp, t_cut, t_null)
  tt <- tt[!is.na(tt)]
  if (length(tt) == 0) NA_real_ else min(tt)
}

# Function: get_first_degradation_method
# Description: Returns a human-readable label for what drove the detected threshold.
#              For v3.2+ verdict tables, reads Consensus_Status and N_Methods_Agree
#              to produce an informative label. Falls back to method-matching for
#              legacy tables.
get_first_degradation_method <- function(vdf) {
  if (is.null(vdf) || nrow(vdf) == 0) return("No Degradation")

  # v3.2+ path
  if ("Consensus_Status" %in% names(vdf)) {
    status <- vdf$Consensus_Status[1]
    n_trig <- if ("N_Methods_Triggered" %in% names(vdf)) vdf$N_Methods_Triggered[1] else NA_integer_
    n_agree <- if ("N_Methods_Agree"    %in% names(vdf)) vdf$N_Methods_Agree[1]     else NA_integer_

    if (is.na(status) || status == "NONE") return("No Degradation")

    if (status == "HARD_FAIL") {
      t_bc  <- if ("Threshold_Observed_BCCeiling" %in% names(vdf)) vdf$Threshold_Observed_BCCeiling[1] else NA_real_
      t_ret <- if ("Threshold_Observed_RetFloor"  %in% names(vdf)) vdf$Threshold_Observed_RetFloor[1]  else NA_real_
      gates <- c(if (!is.na(t_bc))  "BC ceiling" else NULL,
                 if (!is.na(t_ret)) "Retention floor" else NULL)
      return(sprintf("Hard fail — %s", paste(gates, collapse = " + ")))
    }

    if (status == "WARNING_SINGLE") {
      # Identify the specific method that triggered (only one will be non-NA)
      t_cp  <- if ("Threshold_Observed_Changepoint" %in% names(vdf)) vdf$Threshold_Observed_Changepoint[1] else NA_real_
      t_cut <- if ("Threshold_Observed_Cutoff"      %in% names(vdf)) vdf$Threshold_Observed_Cutoff[1]      else NA_real_
      t_null<- if ("Threshold_Observed_NullModel"   %in% names(vdf)) vdf$Threshold_Observed_NullModel[1]   else NA_real_
      t_bc  <- if ("Threshold_Observed_BCCeiling"   %in% names(vdf)) vdf$Threshold_Observed_BCCeiling[1]   else NA_real_
      t_ret <- if ("Threshold_Observed_RetFloor"    %in% names(vdf)) vdf$Threshold_Observed_RetFloor[1]    else NA_real_
      which_fired <- if (!is.na(t_cp))   "Changepoint"      else
                     if (!is.na(t_cut))  "Distance Cutoff"  else
                     if (!is.na(t_null)) "Null Model"        else
                     if (!is.na(t_bc))   "BC Ceiling"        else
                     if (!is.na(t_ret))  "Retention Floor"   else "Unknown"
      n_total <- if (!is.na(n_trig) && n_trig > 0) n_trig else 5L
      return(sprintf("Warning — %s only (1/%d methods)", which_fired, n_total))
    }

    if (status == "CONFIRMED" && !is.na(n_agree)) return(sprintf("Consensus (%d methods agree)", n_agree))
    return("Confirmed")
  }

  # Legacy fallback
  t_cp   <- if ("Threshold_Observed_Changepoint" %in% names(vdf)) vdf$Threshold_Observed_Changepoint[1] else NA_real_
  t_cut  <- if ("Threshold_Observed_Cutoff"      %in% names(vdf)) vdf$Threshold_Observed_Cutoff[1]      else NA_real_
  t_null <- if ("Threshold_Observed_NullModel"   %in% names(vdf)) vdf$Threshold_Observed_NullModel[1]   else NA_real_
  f_thresh <- min(c(t_cp, t_cut, t_null)[!is.na(c(t_cp, t_cut, t_null))])
  if (is.infinite(f_thresh)) return("No Degradation")
  if (!is.na(t_cp)   && f_thresh == t_cp)   return("Changepoint")
  if (!is.na(t_cut)  && f_thresh == t_cut)  return("Distance Cutoff")
  if (!is.na(t_null) && f_thresh == t_null) return("Null Model")
  "Unknown"
}

# ==============================================================================
# SECTION 4: PLOTTING FUNCTIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: plot_dissimilarity_robust (Panel A)
# Description: Plots the Community Dissimilarity (Beta Diversity) vs Trimming.
#              - Red Line: Real Data behavior.
#              - Blue Ribbon: Simulation Baseline (Null Model).
#              - Dynamic Scaling: Zooms in if dissimilarity is low.
# ------------------------------------------------------------------------------
plot_dissimilarity_robust <- function(study_data,
                                      sim_baseline_df,
                                      verdict_df,
                                      max_step_bp,
                                      taxonomic_level,
                                      observed_thresh_bp,
                                      plot_max_step_index = NULL,
                                      warning_thresh_bp = NA_real_) {

  # 1. Determine Increment
  actual_increment <- if (!is.null(study_data$increment)) {
    study_data$increment
  } else {
    unique_trim <- sort(unique(study_data$dissim_data$Trim_BP))
    if (length(unique_trim) > 1) unique_trim[2] - unique_trim[1] else 10
  }

  # 2. Determine Axis Limit (Steps)
  x_axis_max_steps <- if (!is.null(plot_max_step_index)) {
    plot_max_step_index
  } else {
    max_step_bp / actual_increment
  }

  message(paste("Plotting dissimilarity for level", taxonomic_level))

  dd <- study_data$dissim_data
  if (!is_valid_df(dd)) return(ggplot() + labs(title = "A) Dissimilarity", subtitle = "No data") + theme_bw())

  # 3. Filter Real Data
  real_curve <- dd %>%
    dplyr::filter(Level == taxonomic_level, !is.na(Dissimilarity)) %>%
    dplyr::mutate(TrimStep = Trim_BP / actual_increment) %>%
    dplyr::filter(TrimStep <= x_axis_max_steps) %>%
    dplyr::arrange(TrimStep)

  if (nrow(real_curve) == 0) return(ggplot() + labs(title = "A) Dissimilarity", subtitle = "No filtered data") + theme_bw())

  # 5. Plot
  p <- ggplot(real_curve, aes(x = TrimStep, y = Dissimilarity)) +
    geom_line(color = "red", linewidth = 1) +
    theme_bw()

  # 6. Simulation Baseline (WITH IMPUTATION)
  max_y_sim <- 0
  if (is_valid_df(sim_baseline_df)) {
    sim_baseline_filtered <- sim_baseline_df %>%
      dplyr::filter(Level == taxonomic_level) %>%
      dplyr::mutate(TrimStep = Trim_BP / actual_increment) %>%
      dplyr::filter(TrimStep <= x_axis_max_steps)

    if (nrow(sim_baseline_filtered) > 0) {
      
      if ("CIUpper" %in% names(sim_baseline_filtered)) max_y_sim <- max(sim_baseline_filtered$CIUpper, na.rm = TRUE)

      p <- p +
        geom_ribbon(data = sim_baseline_filtered, aes(x = TrimStep, ymin = CILower, ymax = CIUpper),
                    fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
        geom_line(data = sim_baseline_filtered, aes(x = TrimStep, y = MeanDissimilarity),
                  color = "blue", alpha = 0.6, inherit.aes = FALSE)
    }
  }

  # 7. Scales & Labels
  global_max <- max(real_curve$Dissimilarity, max_y_sim, na.rm = TRUE)
  upper_limit <- min(1.0, max(0.1, global_max * 1.1))

  subtitle_text <- "Real Data (Red) vs. Simulated Baseline (Blue)"


  p <- p +
    scale_y_continuous(limits = c(0, upper_limit), expand = expansion(mult = c(0, 0.05)), name = "Bray-Curtis Dissimilarity") +
    scale_x_continuous(limits = c(0, x_axis_max_steps), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title    = paste("A) Community Dissimilarity vs. Trimming (", taxonomic_level, ")", sep = ""),
      subtitle = subtitle_text,
      x        = paste("Trim Step (1 step =", actual_increment, "bp)")
    )

  if (!is.null(verdict_df) && nrow(verdict_df) > 0) {
    tryCatch({
      p <- add_threshold_lines(p, verdict_df$Threshold_Required[1], observed_thresh_bp,
                               actual_increment, warning_thresh_bp = warning_thresh_bp)
    }, error = function(e) {})
  }

  return(p)
}

# ------------------------------------------------------------------------------
# Function: plot_retention_robust (Panel B)
# Description: Plots Taxon Retention (Alpha Diversity) vs Trimming.
#              - Dynamic Scaling: If retention is high (>90%), zooms in the Y-axis.
# ------------------------------------------------------------------------------
plot_retention_robust <- function(study_data,
                                  retention_baseline,
                                  verdict_df,
                                  max_step_bp,
                                  taxonomic_level,
                                  observed_thresh_bp,
                                  plot_max_step_index = NULL,
                                  warning_thresh_bp = NA_real_) {

  # 1. Increment
  actual_increment <- if (!is.null(study_data$increment)) {
    study_data$increment
  } else {
    unique_trim <- sort(unique(study_data$retention_data$Trim_BP))
    if (length(unique_trim) > 1) unique_trim[2] - unique_trim[1] else 10
  }

  # 2. Axis Max
  x_axis_max_steps <- if (!is.null(plot_max_step_index)) {
    plot_max_step_index
  } else {
    max_step_bp / actual_increment
  }

  rd <- study_data$retention_data
  if (!is_valid_df(rd)) return(ggplot() + labs(title = "B) Retention", subtitle = "No data") + theme_bw())

  # 3. Filter Real Data
  real_curve <- rd %>%
    dplyr::filter(Level == taxonomic_level, !is.na(Retention)) %>%
    dplyr::mutate(TrimStep = Trim_BP / actual_increment) %>%
    dplyr::filter(TrimStep <= x_axis_max_steps) %>%
    dplyr::arrange(TrimStep)

  if (nrow(real_curve) == 0) return(ggplot() + labs(title = "B) Retention", subtitle = "No filtered data") + theme_bw())

  # 4. IMPUTE TOTAL LOSS (Drop to 0%)
  last_real_step <- max(real_curve$TrimStep)
  

  # 5. Plot
  p <- ggplot(real_curve, aes(x = TrimStep, y = Retention)) +
    geom_line(color = "red", linewidth = 1) +
    theme_bw()

  # 6. Simulation Baseline (WITH IMPUTATION)
  if (is_valid_df(retention_baseline)) {
    sim_retention_filtered <- retention_baseline %>%
      dplyr::filter(Level == taxonomic_level) %>%
      dplyr::mutate(TrimStep = Trim_BP / actual_increment) %>%
      dplyr::filter(TrimStep <= x_axis_max_steps)

    if (nrow(sim_retention_filtered) > 0) {
      p <- p +
        geom_ribbon(data = sim_retention_filtered, aes(x = TrimStep, ymin = CILower, ymax = CIUpper),
                    fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
        geom_line(data = sim_retention_filtered, aes(x = TrimStep, y = MeanRetention),
                  color = "blue", alpha = 0.6, inherit.aes = FALSE)
    }
  }

  # 7. Dynamic Y-Axis
  min_retention_obs <- min(real_curve$Retention, na.rm = TRUE)
  y_min_limit <- if (min_retention_obs > 90) 80 else 0

  p <- p +
    scale_y_continuous(limits = c(y_min_limit, 101), expand = expansion(mult = c(0, 0)), name = "Taxa Retained (%)") +
    scale_x_continuous(limits = c(0, x_axis_max_steps), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title    = paste("B) Taxon Retention vs. Trimming (", taxonomic_level, ")", sep = ""),
      x        = paste("Trim Step (1 step =", actual_increment, "bp)")
    )

  if (!is.null(verdict_df) && nrow(verdict_df) > 0) {
    tryCatch({
      p <- add_threshold_lines(p, verdict_df$Threshold_Required[1], observed_thresh_bp,
                               actual_increment, warning_thresh_bp = warning_thresh_bp)
    }, error = function(e) {})
  }

  return(p)
}

# ------------------------------------------------------------------------------
# Function: plot_taxon_impact_combined (Panel C)
# Description: Visualizes the specific taxa being affected.
#              - Top: Line chart of the 10 most variable taxa.
#              - Bottom: Bar chart of Top 5 Increased vs Top 5 Decreased taxa.
# ------------------------------------------------------------------------------
plot_taxon_impact_combined <- function(study_data,
                                       verdict_df,
                                       taxonomic_level,
                                       max_step_bp,
                                       observed_thresh_bp,
                                       plot_max_step_index = NULL,
                                       warning_thresh_bp = NA_real_) {

  # 1. Increment & Axis Logic
  actual_increment <- if (!is.null(study_data$increment)) {
    study_data$increment
  } else {
    10 # fallback
  }

  x_axis_max_steps <- if (!is.null(plot_max_step_index)) {
    plot_max_step_index
  } else {
    max_step_bp / actual_increment
  }
  
  max_trim_bp_calculated <- x_axis_max_steps * actual_increment

  message(paste("Creating combined taxon impact visualization for level", taxonomic_level))

  if (is.null(study_data$taxon_impacts)) {
    return(ggplot() + labs(title = "C) Most Impacted Taxa", subtitle = "No taxon impact data") + theme_bw())
  }

  impacted_taxa_raw <- study_data$taxon_impacts$impacted_taxa
  if (!is_valid_df(impacted_taxa_raw)) {
    return(ggplot() + labs(title = "C) Most Impacted Taxa", subtitle = "No impacted taxa data") + theme_bw())
  }

  # Filter Data
  impact_data <- impacted_taxa_raw %>%
    dplyr::filter(Level == taxonomic_level, Trim_BP <= max_trim_bp_calculated) %>%
    dplyr::mutate(TrimStep = Trim_BP / actual_increment)

  # --- PREPARE LINE CHART DATA (With Imputation) ---
  # Find top 10 variable taxa first
  taxa_variability <- impact_data %>%
    dplyr::group_by(Taxon) %>%
    dplyr::summarise(Range = max(Abundance, na.rm = TRUE) - min(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(Range)) %>%
    dplyr::slice_head(n = 10)

  line_chart_data <- impact_data %>%
    dplyr::filter(Taxon %in% taxa_variability$Taxon)

  # IMPUTE: For each taxon in the chart, if data stops early, append rows dropping to 0
  last_real_step <- max(line_chart_data$TrimStep)
  

  # --- PLOT 1: Line Chart ---
  line_chart <- ggplot(line_chart_data, aes(x = TrimStep, y = Abundance, color = Taxon, group = Taxon)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~Taxon, scales = "free_y", ncol = 5) +
    labs(
      title = "Taxa with the largest change in relative abundance",
      x     = paste("Trim Step (1 step =", actual_increment, "bp)"),
      y     = "Mean Relative Abundance"
    ) +
    scale_x_continuous(limits = c(0, x_axis_max_steps), expand = expansion(mult = c(0, 0.02))) +
    theme_bw() +
    theme(legend.position = "none", strip.text = element_text(size = 8), axis.text = element_text(size = 7))

  if (!is.null(verdict_df) && nrow(verdict_df) > 0) {
      tryCatch({
        line_chart <- add_threshold_lines(line_chart, verdict_df$Threshold_Required[1],
                                          observed_thresh_bp, actual_increment,
                                          warning_thresh_bp = warning_thresh_bp)
      }, error = function(e) {})
  }

  # --- PLOT 2: Bar Chart ---
  # For CONFIRMED degradation use confirmed threshold; for WARNING_SINGLE use the warning point;
  # for PASS use max trim (shows full-range abundance change).
  safe_threshold_bp <- if (!is.na(observed_thresh_bp)) {
    observed_thresh_bp
  } else if (!is.na(warning_thresh_bp)) {
    warning_thresh_bp
  } else {
    max_trim_bp_calculated
  }

  # Ensure we don't pick a threshold beyond the plot
  safe_threshold_bp <- min(safe_threshold_bp, max_trim_bp_calculated)

  # Get Initial (0 bp)
  initial_abun <- impact_data %>%
    filter(Trim_BP == 0) %>%
    select(Taxon, InitialAbundance = Abundance)

  # Get Final (at threshold)
  # Logic: If threshold exists in data, use it. If not (degraded), assume 0.
  final_abun_real <- impact_data %>%
    filter(Trim_BP == safe_threshold_bp) %>%
    select(Taxon, FinalAbundance = Abundance)
    
  # Combine
  abundance_changes <- initial_abun %>%
    left_join(final_abun_real, by = "Taxon") %>%
    mutate(
      # If FinalAbundance is NA (meaning the data didn't reach safe_threshold_bp),
      # replace with 0.0 assuming total loss.
      FinalAbundance = ifelse(is.na(FinalAbundance), 0.0, FinalAbundance),
      AbundanceChange = FinalAbundance - InitialAbundance
    )

  top_increased <- abundance_changes %>% arrange(desc(AbundanceChange)) %>% slice_head(n = 5) %>% mutate(ChangeType = "Increased")
  top_decreased <- abundance_changes %>% arrange(AbundanceChange) %>% slice_head(n = 5) %>% mutate(ChangeType = "Decreased")
  bar_chart_data <- bind_rows(top_increased, top_decreased)

  if (nrow(bar_chart_data) > 0) {
    bar_chart <- ggplot(bar_chart_data, aes(x = reorder(Taxon, AbundanceChange), y = AbundanceChange)) +
      geom_col(aes(fill = AbundanceChange > 0), width = 0.7) +
      scale_fill_manual(values = c("FALSE" = "darkred", "TRUE" = "darkgreen"), labels = c("Decrease", "Increase"), name = "Direction") +
      coord_flip() +
      labs(title = paste("Top 5 Increased + Top 5 Decreased", taxonomic_level, "s"),
           subtitle = paste("Change from 0 to", safe_threshold_bp, "bp"),
           x = taxonomic_level, y = "Change in Relative Abundance") +
      theme_bw() + theme(legend.position = "bottom")
  } else {
    bar_chart <- ggplot() + labs(title = "Top Changed Taxa", subtitle = "No abundance change data") + theme_bw()
  }

  combined_plot <- wrap_plots(line_chart, bar_chart, ncol = 1, heights = c(2, 1)) +
    plot_annotation(title = paste("C) Most Impacted", taxonomic_level, "s"))

  return(combined_plot)
}

# ------------------------------------------------------------------------------
# Function: plot_core_taxa_robust (Panel D)
# Description: Plots the Core Taxa (High Prevalence, High Abundance).
#              - Fallback: If no core taxa exist, shows the practical study inclusion guide.
# ------------------------------------------------------------------------------
plot_core_taxa_robust <- function(study_data,
                                  verdict_df,
                                  max_step_bp,
                                  taxonomic_level,
                                  observed_thresh_bp,
                                  plot_max_step_index = NULL,
                                  warning_thresh_bp = NA_real_) {

  # 1. Increment & Axis Logic
  actual_increment <- if (!is.null(study_data$increment)) {
    study_data$increment
  } else {
    10 # fallback
  }

  x_axis_max_steps <- if (!is.null(plot_max_step_index)) {
    plot_max_step_index
  } else {
    max_step_bp / actual_increment
  }
  
  max_trim_bp_calculated <- x_axis_max_steps * actual_increment

  message(paste("Creating core taxa plot for level:", taxonomic_level))

  # Check if taxon_impacts exists at all
  if (is.null(study_data$taxon_impacts)) {
    message("   -> No taxon impacts data - showing practical guide instead")
    return(create_study_specific_practical_guide(study_data, verdict_df))
  }

  # Config
  n_taxa_to_plot <- if (exists("TOP_N_CORE_TAXA")) TOP_N_CORE_TAXA else 5

  # ONLY try core_taxa - no fallback to impacted_taxa
  core_data <- NULL
  if (!is.null(study_data$taxon_impacts$core_taxa) && nrow(study_data$taxon_impacts$core_taxa) > 0) {
     core_data <- study_data$taxon_impacts$core_taxa %>% 
       filter(Level == taxonomic_level, Trim_BP <= max_trim_bp_calculated)
     
     # Pick Top N by mean abundance
     if (nrow(core_data) > 0) {
       top_cores <- core_data %>% 
         group_by(Taxon) %>% 
         summarise(M = mean(Abundance, na.rm = TRUE), .groups = "drop") %>% 
         arrange(desc(M)) %>% 
         slice_head(n = n_taxa_to_plot)
       core_data <- core_data %>% filter(Taxon %in% top_cores$Taxon)
     }
  }

  # If no core taxa, show practical guide
  if (is.null(core_data) || nrow(core_data) == 0) {
    message("   -> No core taxa available - showing practical guide instead")
    return(create_study_specific_practical_guide(study_data, verdict_df))
  }

  # Add Steps
  core_data <- core_data %>% mutate(TrimStep = Trim_BP / actual_increment)

  # IMPUTE: Drop to 0 if data ends early
  last_real_step <- max(core_data$TrimStep)
  

  # Plot
  p <- ggplot(core_data, aes(x = TrimStep, y = Abundance, color = Taxon, group = Taxon)) +
    geom_line(linewidth = 1) +
    scale_x_continuous(limits = c(0, x_axis_max_steps), expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~Taxon, scales = "free_y", ncol = 3) +
    labs(
      title = paste("D) Core", taxonomic_level, "s"),
      subtitle = paste0("Top ", n_taxa_to_plot, " abundant taxa (prevalence-filtered)"),
      x = paste("Trim Step (1 step =", actual_increment, "bp)"),
      y = "Mean Relative Abundance"
    ) +
    theme_bw() + theme(legend.position = "none", strip.text = element_text(size = 8))

  if (!is.null(verdict_df) && nrow(verdict_df) > 0) {
    tryCatch({
      p <- add_threshold_lines(p, verdict_df$Threshold_Required[1], observed_thresh_bp,
                               actual_increment, warning_thresh_bp = warning_thresh_bp)
    }, error = function(e) {})
  }

  return(p)
}

# ------------------------------------------------------------------------------
# Helper Function: create_study_specific_practical_guide
# Description: Creates a simplified practical guide for a single study (Panel D fallback)
# ------------------------------------------------------------------------------
create_study_specific_practical_guide <- function(study_data, verdict_df) {
  
  if (is.null(verdict_df) || nrow(verdict_df) == 0) {
    return(ggplot() + 
      labs(title = "D) Core Taxa", subtitle = "No data available") + 
      theme_bw())
  }
  
  study_name <- if (!is.null(study_data$study_name)) study_data$study_name else "Unknown"
  primer <- if (!is.null(study_data$primer_name)) study_data$primer_name else "Unknown"
  target_trim <- if (!is.null(study_data$target_trim_bp)) study_data$target_trim_bp else NA
  
  # Define dissimilarity warning threshold (configurable)
  DISSIM_WARNING_THRESHOLD <- if (exists("DISSIM_WARNING_THRESHOLD")) {
    DISSIM_WARNING_THRESHOLD 
  } else {
    0.08  # Conservative default
  }
  
  # Analyze each level
  level_summary <- verdict_df %>%
    mutate(
      First_Threshold = pmap_dbl(
        list(
          Threshold_Observed_Changepoint, Threshold_Observed_Cutoff, Threshold_Observed_NullModel,
          if ("Confirmed_Threshold_BP" %in% names(.)) Confirmed_Threshold_BP else rep(NA_real_, n()),
          if ("Consensus_Status"       %in% names(.)) Consensus_Status       else rep(NA_character_, n())
        ),
        compute_first_threshold
      ),
      Status = case_when(
        is.na(First_Threshold) &
          "Consensus_Status" %in% names(.) &
          Consensus_Status == "WARNING_SINGLE"              ~ "CAUTION",
        is.na(First_Threshold)                              ~ "PASS",
        First_Threshold >= Threshold_Required               ~ "PASS",
        First_Threshold >= (Threshold_Required * 0.85)      ~ "MARGINAL",
        TRUE                                                ~ "FAIL"
      ),
      StatusSymbol = case_when(
        Status == "PASS"    ~ "✓",
        Status == "CAUTION" ~ "~",
        Status == "MARGINAL" ~ "⚠",
        TRUE                ~ "✗"
      ),
      Gap = ifelse(is.na(First_Threshold), NA, First_Threshold - Threshold_Required)
    ) %>%
    select(Level, Status, StatusSymbol, First_Threshold, Threshold_Required, Gap)
  
  # Check dissimilarity at target trim point for PASS levels
  dissim_warnings <- c()
  
  if (!is.null(study_data$dissim_data) && !is.na(target_trim) && target_trim > 0) {
    for (i in 1:nrow(level_summary)) {
      if (level_summary$Status[i] == "PASS") {
        level_name <- level_summary$Level[i]
        
        # Get dissimilarity at target trim point
        dissim_at_target <- study_data$dissim_data %>%
          filter(Level == level_name, Trim_BP == target_trim) %>%
          pull(Dissimilarity)
        
        if (length(dissim_at_target) > 0 && !is.na(dissim_at_target[1])) {
          dissim_value <- dissim_at_target[1]
          
          if (dissim_value >= DISSIM_WARNING_THRESHOLD) {
            warning_text <- sprintf("%s: β-div = %.3f", level_name, dissim_value)
            dissim_warnings <- c(dissim_warnings, warning_text)
            
            # Update status symbol to warning
            level_summary$StatusSymbol[i] <- "⚠"
            level_summary$Status[i] <- "PASS*"
          }
        }
      }
    }
  }
  
  # Overall verdict logic
  n_total <- nrow(level_summary)
  n_pass <- sum(grepl("PASS", level_summary$Status), na.rm = TRUE)  # Includes PASS*
  n_marginal <- sum(level_summary$Status == "MARGINAL", na.rm = TRUE)
  n_fail <- sum(level_summary$Status == "FAIL", na.rm = TRUE)
  
  pass_rate <- n_pass / n_total
  moderate_rate <- (n_pass + n_marginal) / n_total
  
  # Adjust overall verdict if dissimilarity warnings present
  has_dissim_warnings <- length(dissim_warnings) > 0
  
  overall_verdict <- if (n_total == 1) {
    if (level_summary$Status[1] == "PASS" && !has_dissim_warnings) {
      "RECOMMENDED (single level)"
    } else if (level_summary$Status[1] == "PASS*" || level_summary$Status[1] == "MARGINAL") {
      "MARGINAL (review β-diversity shift)"
    } else {
      "NOT RECOMMENDED"
    }
  } else if (has_dissim_warnings && pass_rate < 1.0) {
    "MARGINAL - Gradual degradation detected"
  } else if (pass_rate >= 0.6) {
    if (has_dissim_warnings) {
      "RECOMMENDED - but review β-diversity"
    } else {
      "RECOMMENDED"
    }
  } else if (moderate_rate >= 0.6) {
    "MARGINAL - Review carefully"
  } else {
    "NOT RECOMMENDED"
  }
  
  # Format table with better alignment
  table_text <- level_summary %>%
    mutate(
      Display = sprintf("%-10s %s %-10s   %s", 
                       Level, 
                       StatusSymbol,
                       Status,
                       ifelse(is.na(First_Threshold), "No degr.", 
                              sprintf("Degr: %4.0f bp", First_Threshold))
      )
    ) %>%
    pull(Display) %>%
    paste(collapse = "\n")
  
  # Build warning section if needed
  warning_section <- if (has_dissim_warnings) {
    paste(
      "",
      "  ⚠ BETA DIVERSITY WARNINGS:",
      paste0("    ", dissim_warnings, collapse = "\n"),
      sprintf("    (Threshold: β-div > %.2f at target trim)", DISSIM_WARNING_THRESHOLD),
      "",
      sep = "\n"
    )
  } else {
    ""
  }
  
  guide_text <- paste(
    "═══════════════════════════════════════════════════════",
    sprintf("  STUDY: %s", study_name),
    sprintf("  PRIMER: %s", primer),
    "═══════════════════════════════════════════════════════",
    "",
    "  QUALITY ASSESSMENT BY TAXONOMIC LEVEL",
    "",
    table_text,
    warning_section,
    "───────────────────────────────────────────────────────",
    sprintf("  OVERALL: %s", overall_verdict),
    "───────────────────────────────────────────────────────",
    "",
    sprintf("  Pass: %d | Marginal: %d | Fail: %d  (of %d levels)", 
            n_pass, n_marginal, n_fail, n_total),
    "",
    if (!is.na(target_trim) && target_trim > 0) {
      sprintf("  Target trim: %d bp", target_trim)
    } else {
      ""
    },
    "",
    "  Note: Core taxa data unavailable for this study.",
    "  This panel shows quality assessment instead.",
    if (has_dissim_warnings) {
      "  * = Passed but exceeds β-diversity threshold"
    } else {
      ""
    },
    "",
    "═══════════════════════════════════════════════════════",
    sep = "\n"
  )
  
  # Match the style of other panels
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
            label = guide_text,
            hjust = 0.5, vjust = 0.5, 
            size = 3.2, family = "mono", lineheight = 1.05) +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      title = "D) Study Quality Assessment",
      subtitle = "Core taxa unavailable - showing pass/fail status"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.0, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.0, size = 10),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

# ==============================================================================
# SECTION 5: IMPROVED MASTER SUMMARY FUNCTIONS
# These functions generate the visualizations for the "Master Summary Report".
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: create_trimming_recommendations
# Description: Generates a Forest Plot summarizing the "Safe Trimming Thresholds"
#              across all studies for each taxonomic level.
#              - Red Dot: Median Threshold.
#              - Blue Bar: Interquartile Range (IQR).
#              - Text: Pass Rate %.
# ------------------------------------------------------------------------------
create_trimming_recommendations <- function() {
    message("--- Creating trimming recommendations ---")

    if (nrow(master_verdicts) == 0) {
        message("⚠️ No verdict data available for recommendations")
        return(NULL)
    }

    # Summarize verdicts by Level
    consensus_data <- master_verdicts %>%
        mutate(
            # Use confirmed threshold from consensus voting (v3.2+), or legacy min()
            First_Threshold = pmap_dbl(
              list(
                Threshold_Observed_Changepoint, Threshold_Observed_Cutoff, Threshold_Observed_NullModel,
                if ("Confirmed_Threshold_BP" %in% names(.)) Confirmed_Threshold_BP else rep(NA_real_, n()),
                if ("Consensus_Status"       %in% names(.)) Consensus_Status       else rep(NA_character_, n())
              ),
              compute_first_threshold
            ),
            # PASS if no confirmed degradation OR confirmed degradation is beyond required trim.
            # WARNING_SINGLE (First_Threshold = NA from compute_first_threshold) counts as PASS here
            # since the lone trigger was not corroborated — it is flagged separately in the subtitle.
            Pass_Status = ifelse(is.na(First_Threshold) | First_Threshold >= Threshold_Required, 1, 0)
        ) %>%
        group_by(Level) %>%
        summarise(
            Total_Studies = n(),
            Pass_Count = sum(Pass_Status, na.rm = TRUE),
            Pass_Rate = round(Pass_Count / Total_Studies * 100, 1),
            Q25_Threshold = quantile(First_Threshold, 0.25, na.rm = TRUE),
            Median_Threshold = median(First_Threshold, na.rm = TRUE),
            Q75_Threshold = quantile(First_Threshold, 0.75, na.rm = TRUE),
            Conservative_Threshold = Q25_Threshold,
            .groups = "drop"
        ) %>%
        arrange(desc(Pass_Rate), desc(Median_Threshold)) %>%
        mutate(
            Recommendation = case_when(
                Pass_Rate >= 85 & Median_Threshold >= 150 ~ "HIGHLY RECOMMENDED",
                Pass_Rate >= 75 & Median_Threshold >= 100 ~ "RECOMMENDED",
                Pass_Rate >= 60 ~ "USE WITH CAUTION",
                TRUE ~ "NOT RECOMMENDED"
            ),
            Safe_Trim_Length = pmax(Conservative_Threshold - 20, 50)
        )

    # Generate Forest Plot
    p <- consensus_data %>%
        ggplot(aes(x = reorder(Level, Median_Threshold))) +
        geom_errorbar(aes(ymin = Q25_Threshold, ymax = Q75_Threshold),
                     width = 0.3, linewidth = 2, color = "steelblue") +
        geom_point(aes(y = Median_Threshold), size = 4, color = "darkred") +
        geom_text(aes(y = Median_Threshold + 20,
                     label = paste0(round(Median_Threshold), " bp")),
                 size = 3, fontface = "bold") +
        geom_text(aes(y = 20,
                     label = paste0(Pass_Rate, "% pass rate")),
                 size = 3, color = "darkgreen") +

        coord_flip() +
        labs(
            title = "Safe Trimming Recommendations by Taxonomic Level",
            subtitle = "Red dot = median threshold, blue bar = IQR (25th-75th percentile)",
            x = "Taxonomic Level",
            y = "Safe Trimming Threshold (bp)"
        ) +

        theme_bw() +
        theme(
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11)
        )

    return(p)
}

# ------------------------------------------------------------------------------
# Function: create_amplicon_alignment_plot
# Description: Visualizes the alignment of all studies relative to the 16S gene.
#              - Handles both "Study Mode" (global alignment) and "Primer Mode" (local).
#              - Automatically detects Coordinate System (E. coli vs SILVA).
#              - Now reads consensus from consensusregioninfo.csv to avoid recalculation.
# ------------------------------------------------------------------------------
create_amplicon_alignment_plot <- function(master_manifest, master_verdicts, taxonomic_level = "Family") {
    message(paste("--- Creating amplicon alignment visualization for:", taxonomic_level, "---"))

    if (nrow(master_verdicts) == 0 || nrow(master_manifest) == 0) {
        message("[WARN] Master verdict or manifest data is missing. Cannot generate plot.")
        return(NULL)
    }

    overlaps_consensus <- function(seg_start, seg_end, cons_start, cons_end) {
        return(seg_start < cons_end && seg_end > cons_start)
    }

    # --- 1. DYNAMIC COORDINATE LOADING ---
    study_coords_path <- "output/intermediate/study_alignment_coords.csv"
    primer_coords_path <- "output/intermediate/primer_coords_phase1_output.csv"

    target_path <- NULL
    coords_source <- ""

    if (file.exists(study_coords_path)) {
        target_path <- study_coords_path
        coords_source <- "study"
        message("   -> Found Study Mode coordinates: ", target_path)
    } else if (file.exists(primer_coords_path)) {
        target_path <- primer_coords_path
        coords_source <- "primer"
        message("   -> Found Primer Mode coordinates: ", target_path)
    } else {
        message("[WARN] No alignment coordinates file found at standard locations.")
        return(NULL)
    }

    # Read coords file with error handling
    primer_coords <- tryCatch({
        pc <- readr::read_csv(target_path, show_col_types = FALSE)

        # Normalize column names if coming from primer tool
        if (coords_source == "primer") {
            if (!("ref_start" %in% names(pc)) && "primer_start" %in% names(pc)) {
                pc <- pc %>% dplyr::rename(ref_start = primer_start, ref_end = primer_end)
            }
        }
        message(paste("   -> Loaded", nrow(pc), "coordinate records"))
        pc
    }, error = function(e) {
        message("[WARN] Error reading coords file: ", conditionMessage(e))
        return(NULL)
    })

    if (is.null(primer_coords) || nrow(primer_coords) == 0) {
        message("[WARN] No primer coordinates found or file is empty.")
        return(NULL)
    }

    # --- 2. PREPARE AMPLICON DATA ---
    amplicon_data <- master_verdicts %>%
        dplyr::filter(Level == taxonomic_level) %>%
        {
            if (coords_source == "study") {
                dplyr::left_join(., primer_coords, by = c("Study" = "study_name"))
            } else {
                dplyr::left_join(., primer_coords, by = c("Primer" = "primer_name"))
            }
        } %>%
        dplyr::mutate(
            # For the alignment map we need the raw trigger position regardless of consensus
            # status — even a WARNING_SINGLE should split the segment visually.
            # Use confirmed threshold if available; fall back to raw min of individual methods.
            # This is intentionally different from compute_first_threshold(), which returns NA
            # for WARNING_SINGLE to suppress the red line on individual study plots.
            # Inner threshold (most conservative = lowest bp = red band boundary)
            First_Threshold = purrr::pmap_dbl(
                list(
                  if ("Threshold_Inner_BP"    %in% names(.)) Threshold_Inner_BP    else rep(NA_real_, dplyr::n()),
                  if ("Confirmed_Threshold_BP" %in% names(.)) Confirmed_Threshold_BP else rep(NA_real_, dplyr::n()),
                  if ("Consensus_Status"       %in% names(.)) Consensus_Status       else rep(NA_character_, dplyr::n()),
                  Threshold_Observed_Changepoint, Threshold_Observed_Cutoff, Threshold_Observed_NullModel
                ),
                function(inner_bp, confirmed_bp, consensus_status, cp, cut, null) {
                  if (!is.na(consensus_status) && consensus_status %in% c("HARD_FAIL", "CONFIRMED"))
                    return(if (!is.na(inner_bp)) inner_bp else confirmed_bp)
                  raw_vals <- c(cp, cut, null)
                  raw_vals <- raw_vals[!is.na(raw_vals)]
                  if (length(raw_vals) == 0) NA_real_ else min(raw_vals)
                }
            ),
            # Outer threshold (most lenient = highest bp = caution band boundary)
            Outer_Threshold = purrr::pmap_dbl(
                list(
                  if ("Threshold_Outer_BP"    %in% names(.)) Threshold_Outer_BP    else rep(NA_real_, dplyr::n()),
                  if ("Consensus_Status"       %in% names(.)) Consensus_Status       else rep(NA_character_, dplyr::n()),
                  Threshold_Observed_Changepoint, Threshold_Observed_Cutoff, Threshold_Observed_NullModel
                ),
                function(outer_bp, consensus_status, cp, cut, null) {
                  if (!is.na(consensus_status) && consensus_status %in% c("HARD_FAIL", "CONFIRMED"))
                    return(outer_bp)
                  # For WARNING_SINGLE: outer = inner (single band, no two-stage)
                  raw_vals <- c(cp, cut, null)
                  raw_vals <- raw_vals[!is.na(raw_vals)]
                  if (length(raw_vals) == 0) NA_real_ else max(raw_vals)
                }
            ),
Is_Warning_Only = !is.na(Consensus_Status) & Consensus_Status == "WARNING_SINGLE",
            Amplicon_Start = ref_start,
            Amplicon_End = ref_end,
            Head_Proportion = if ("Head_Proportion" %in% names(.)) pmax(0, pmin(1, Head_Proportion)) else 0.5,
            Study_Label = dplyr::case_when(
                coords_source == "study" ~ as.character(Study),
                TRUE ~ as.character(Primer)
            )
        ) %>%
        dplyr::arrange(Study) %>%
        dplyr::mutate(Study_Position = row_number(),
                      Label_X        = -(max(ref_end, na.rm = TRUE) * 1.05 * 0.02) * 2.0) %>%
        dplyr::ungroup() %>%
        dplyr::filter(!is.na(Amplicon_Start) & !is.na(Amplicon_End))

    if (nrow(amplicon_data) == 0) {
        message(paste("[WARN] No data with valid coordinates found for level:", taxonomic_level))
        return(NULL)
    }

    message(paste("   -> Plotting", nrow(amplicon_data), "studies"))

    # --- 3. CONSENSUS REGION LOADING FROM FILE (NEW LOGIC) ---
    consensus_file <- file.path(OUTDIR, "intermediate/consensusregioninfo.csv")
    outlier_studies <- character(0)
    consensus_start <- NA
    consensus_end <- NA
    consensus_source <- "Not Found"
    
    if (file.exists(consensus_file)) {
        consensus_info <- suppressMessages(readr::read_csv(consensus_file, show_col_types = FALSE))
        
        if (nrow(consensus_info) > 0) {
            consensus_start <- consensus_info$ConsensusStart[1]
            consensus_end <- consensus_info$ConsensusEnd[1]
            consensus_source <- "Consensus File"
            
            if (!is.na(consensus_info$OutlierStudies[1]) && consensus_info$OutlierStudies[1] != "") {
                outlier_studies <- unlist(strsplit(consensus_info$OutlierStudies[1], ";"))
                outlier_studies <- trimws(outlier_studies)
                message(paste("   -> Loaded consensus from file:", consensus_start, "-", consensus_end))
                message(paste("   -> Outliers to exclude:", paste(outlier_studies, collapse = ", ")))
            }
        }
    }
    
    # Fallback: Recalculate if file doesn't exist or is invalid
    if (is.na(consensus_start) || is.na(consensus_end) || consensus_start >= consensus_end) {
        message("   [WARN] No valid consensus file found. Recalculating from scratch...")
        source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_consensus.R"))
        
        consres <- find_largest_overlapping_clique(
            starts = amplicon_data$Amplicon_Start,
            ends = amplicon_data$Amplicon_End,
            study_names = amplicon_data$Study_Label,
            min_overlap = 50
        )
        
        consensus_start <- consres$start
        consensus_end <- consres$end
        outlier_studies <- consres$excluded_studies
        consensus_source <- "Intersection (Calculated)"
    }

    message(paste("   -> Using Consensus:", consensus_start, "-", consensus_end, "(", consensus_source, ")"))

    # --- 4. MARK OUTLIERS IN AMPLICON DATA ---
    amplicon_data <- amplicon_data %>%
        dplyr::mutate(
            Is_Outlier = Study_Label %in% outlier_studies | (if("Is_Outlier" %in% names(.)) Is_Outlier else FALSE)
        )

# --- 5. CREATE PLOT SEGMENTS ---
    add_seg <- function(segs, study, pos, start, end, type) {
        if (is.na(start) || is.na(end) || end <= start) return(segs)
        dplyr::bind_rows(segs, tibble::tibble(
            Study = study, Study_Position = pos,
            Segment_Start = start, Segment_End = end, Segment_Type = type))
    }

    COL_OUT   <- "Outside consensus"
    COL_GREEN <- "Safe - within consensus"
    COL_WARN  <- "Caution (1-method warning)"
    COL_RED   <- "Degraded (confirmed)"
    COL_GREY  <- "Outlier (Excluded)"

    plot_segments <- tibble::tibble()

    for (i in seq_len(nrow(amplicon_data))) {
        row             <- amplicon_data[i, ]
        is_outlier      <- isTRUE(row$Is_Outlier)
        is_warning_only <- isTRUE(row$Is_Warning_Only)
        amp_s     <- row$Amplicon_Start
        amp_e     <- row$Amplicon_End
        lbl       <- row$Study_Label
        pos       <- row$Study_Position
        cs        <- consensus_start
        ce        <- consensus_end
        dg_raw    <- row$First_Threshold
        outer_raw <- if (!is.na(row$Outer_Threshold)) row$Outer_Threshold else dg_raw
        two_stage <- !is.na(outer_raw) && !is.na(dg_raw) && (outer_raw - dg_raw) > 50

        # OUTLIER
        if (is_outlier) {
            plot_segments <- add_seg(plot_segments, lbl, pos, amp_s, amp_e, COL_GREY)
            next
        }

        # NO DEGRADATION — green inside consensus, light orange outside
        if (is.na(dg_raw) || dg_raw <= 0) {
            lbl_star <- paste0(lbl, " *")
            amplicon_data$Study_Label[i] <- lbl_star
            plot_segments <- add_seg(plot_segments, lbl_star, pos, amp_s,          min(amp_e, cs), COL_OUT)
            plot_segments <- add_seg(plot_segments, lbl_star, pos, max(amp_s, cs), min(amp_e, ce), COL_GREEN)
            plot_segments <- add_seg(plot_segments, lbl_star, pos, max(amp_s, ce), amp_e,          COL_OUT)
            next
        }

        # Compute absolute reference positions of degradation boundaries
        head_trim_inner <- round(dg_raw    * row$Head_Proportion)
        tail_trim_inner <- dg_raw    - head_trim_inner
        head_trim_outer <- round(outer_raw * row$Head_Proportion)
        tail_trim_outer <- outer_raw - head_trim_outer

        # Absolute positions: NA means that end contributes no trim
        dg_L_inner <- if (head_trim_inner > 0) amp_s + head_trim_inner else NA_real_
        dg_R_inner <- if (tail_trim_inner > 0) amp_e - tail_trim_inner else NA_real_
        dg_L_outer <- if (head_trim_outer > 0) amp_s + head_trim_outer else NA_real_
        dg_R_outer <- if (tail_trim_outer > 0) amp_e - tail_trim_outer else NA_real_

        # Determine where degradation falls relative to consensus
        L_in_cons <- !is.na(dg_L_inner) && dg_L_inner > cs && dg_L_inner < ce
        R_in_cons <- !is.na(dg_R_inner) && dg_R_inner > cs && dg_R_inner < ce

        # Choose degradation colour
        col_dg <- if (is_warning_only) COL_WARN else COL_RED

        if (!L_in_cons && !R_in_cons) {
            # ── DEGRADATION ENTIRELY IN OVERHANG(S) ──────────────────────
            # Degradation detected before consensus — consensus itself cannot
            # be considered safe, so paint it with the degradation colour too.

            # Left overhang
            if (!is.na(dg_L_inner) && dg_L_inner <= cs) {
                # Degradation starts in left overhang
                plot_segments <- add_seg(plot_segments, lbl, pos, amp_s,       dg_L_inner,     COL_OUT)
                plot_segments <- add_seg(plot_segments, lbl, pos, dg_L_inner,  min(amp_e, cs), col_dg)
            } else {
                plot_segments <- add_seg(plot_segments, lbl, pos, amp_s,       min(amp_e, cs), COL_OUT)
            }

            # Consensus — degraded (signal already present in overhang)
            plot_segments <- add_seg(plot_segments, lbl, pos, max(amp_s, cs), min(amp_e, ce), col_dg)

            # Right overhang
            if (!is.na(dg_R_inner) && dg_R_inner >= ce) {
                # Degradation continues into right overhang
                plot_segments <- add_seg(plot_segments, lbl, pos, max(amp_s, ce), dg_R_inner, col_dg)
                plot_segments <- add_seg(plot_segments, lbl, pos, dg_R_inner,     amp_e,      COL_OUT)
            } else {
                plot_segments <- add_seg(plot_segments, lbl, pos, max(amp_s, ce), amp_e,      COL_OUT)
            }

        } else {
            # ── DEGRADATION INSIDE CONSENSUS ─────────────────────────────

            # Left overhang always COL_OUT
            plot_segments <- add_seg(plot_segments, lbl, pos, amp_s, min(amp_e, cs), COL_OUT)

            # Green from cs to where degradation starts (use outer boundary for two-stage)
            safe_end_L <- if (L_in_cons) {
                if (two_stage && !is.na(dg_L_outer)) dg_L_outer else dg_L_inner
            } else cs
            plot_segments <- add_seg(plot_segments, lbl, pos,
                                     max(amp_s, cs), min(amp_e, safe_end_L), COL_GREEN)

            # Two-stage caution band (outer → inner on left)
            if (two_stage && L_in_cons && !is.na(dg_L_outer) && !is.na(dg_L_inner)) {
                plot_segments <- add_seg(plot_segments, lbl, pos,
                                         max(amp_s, dg_L_outer), min(amp_e, dg_L_inner), COL_WARN)
            }

            # Core degraded zone
            core_L <- if (L_in_cons && !is.na(dg_L_inner)) dg_L_inner else max(amp_s, cs)
            core_R <- if (R_in_cons && !is.na(dg_R_inner)) dg_R_inner else min(amp_e, ce)
            plot_segments <- add_seg(plot_segments, lbl, pos, core_L, core_R, col_dg)

            # Two-stage caution band (inner → outer on right)
            if (two_stage && R_in_cons && !is.na(dg_R_inner) && !is.na(dg_R_outer)) {
                plot_segments <- add_seg(plot_segments, lbl, pos,
                                         max(amp_s, dg_R_inner), min(amp_e, dg_R_outer), COL_WARN)
            }

            # Green from degradation end to ce (use outer boundary for two-stage)
            safe_start_R <- if (R_in_cons) {
                if (two_stage && !is.na(dg_R_outer)) dg_R_outer else dg_R_inner
            } else ce
            plot_segments <- add_seg(plot_segments, lbl, pos,
                                     max(amp_s, safe_start_R), min(amp_e, ce), COL_GREEN)

            # Right overhang always COL_OUT
            plot_segments <- add_seg(plot_segments, lbl, pos, max(amp_s, ce), amp_e, COL_OUT)
        }
    }

            # --- 6. DEFINE V-REGIONS (Coordinate System Detection) ---
    # V-region positions in SILVA 138.1 SSU alignment column coordinates
    # (i.e. which(chars != "-") space, matching the ref_start/ref_end values in
    # study_alignment_coords.csv produced by secat_mapping.R via DECIPHER::AlignProfiles)
    # E. coli K-12 16S mapped through SILVA 138.1 full-length alignment.
    # Source: Yarza et al. 2014 (Nat Rev Microbiol) positions projected onto SILVA 138.1.
    # These are approximate; the consensus region boundaries in the data (~13902-23446)
    # confirm V4-V5 falls in this range, consistent with the values below.
    v_regions_silva_aln <- data.frame(
        Region = c("V1",   "V2",   "V3",   "V4",   "V5",   "V6",   "V7",   "V8",   "V9"),
        Start  = c(1044,   2315,   6428,  13902,  16073,  19819,  23820,  27469,  28954),
        End    = c(2315,   4343,   9906,  16073,  19819,  23820,  27469,  28954,  43180)
    )

    # E. coli bp coordinates (only used as fallback if data_max < 5000)
    v_regions_ecoli <- data.frame(
        Region = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
        Start  = c(69,   137,  338,  576,  822,  986,  1117, 1253, 1435),
        End    = c(99,   242,  533,  681,  879,  1043, 1173, 1313, 1465)
    )

    data_max <- max(amplicon_data$Amplicon_End, na.rm = TRUE)

    if (data_max > 5000) {
        v_regions   <- v_regions_silva_aln
        axis_label  <- "SILVA 138.1 SSU Alignment Column"
    } else {
        v_regions   <- v_regions_ecoli
        axis_label  <- "E. coli 16S Position (bp)"
    }

    plot_min <- min(amplicon_data$Amplicon_Start, na.rm=TRUE)
    v_regions <- v_regions %>%
      dplyr::filter(End >= plot_min & Start <= data_max)

    n_studies <- nrow(amplicon_data)

    # --- 7. DYNAMIC SCALING FOR LABELS ---
    max_x_val <- max(amplicon_data$Amplicon_End, na.rm=TRUE)
    max_x_limit <- ceiling((max_x_val * 1.05) / 100) * 100
    label_x_pos <- -(max_x_limit * 0.02)
    left_limit <- label_x_pos * 2.0

    # --- 8. BUILD PLOT ---
    p <- ggplot2::ggplot() +
        # Consensus Box
        ggplot2::annotate("rect", xmin = consensus_start, xmax = consensus_end, ymin = 0, ymax = n_studies + 1, fill = "#e0f3ff", alpha = 0.7) +
        ggplot2::annotate("text", x = (consensus_start + consensus_end) / 2, y = n_studies + 0.7, label = "Consensus Region", size = 3.5, fontface = "bold", color = "#005a9e") +
        # V-Regions
        ggplot2::geom_text(data = v_regions, aes(x = (Start + End) / 2, y = n_studies + 1.2, label = Region), size = 3, fontface = "bold", color = "gray50") +
        ggplot2::geom_vline(data = v_regions, aes(xintercept = Start), linetype = "dashed", alpha = 0.5, color = "gray70") +
        # Study Segments
        ggplot2::geom_segment(data = plot_segments, aes(x = Segment_Start, xend = Segment_End, y = Study_Position, yend = Study_Position, color = Segment_Type), linewidth = 6, alpha = 1.0) +
        # Study Labels
        ggplot2::geom_text(data = amplicon_data, aes(x = Label_X, y = Study_Position, label = Study_Label), hjust = 1, size = 3, family = "sans") +
        ggplot2::scale_color_manual(
            name = "Segment Type",
            values = c(
                "Safe - within consensus"    = "#27ae60",   # green
                "Outside consensus"          = "#f5cba7",   # light orange
                "Caution (1-method warning)" = "#e67e22",   # dark orange
                "Degraded (confirmed)"       = "#e74c3c",   # red
                "Outlier (Excluded)"         = "#95a5a6"    # grey
            ),
            guide = guide_legend(override.aes = list(linewidth = 4))
        ) +
        ggplot2::scale_x_continuous(
            breaks = scales::pretty_breaks(n = 10),
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(limits = c(0, n_studies + 1.5), expand = c(0, 0)) +
        ggplot2::labs(
            title    = "Overall Amplicon Alignment and Safe Trimming Regions",
            subtitle = paste("Analysis at the", taxonomic_level, "Level -", n_studies, "studies"),
            x        = axis_label,
            y        = NULL,
            caption  = "* No degradation detected up to tested trim points"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            legend.position = "top", legend.title = element_text(face = "bold"),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.title.x = element_text(size = 12, face = "bold"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 100, unit = "pt")
        ) +
        ggplot2::coord_cartesian(clip = "off", xlim = c(0, max_x_limit), expand = FALSE)

    return(p)
}

# ------------------------------------------------------------------------------
# Function: create_practical_guide
# Description: Generates a practical study inclusion guide showing which studies
#              can be safely included at each taxonomic level.
# ------------------------------------------------------------------------------
create_practical_guide <- function() {
    message("--- Creating practical study inclusion guide ---")
    
    if (nrow(master_verdicts) == 0) {
        message("⚠️ No verdict data available for practical guide")
        return(NULL)
    }
    
    # Define dissimilarity warning threshold (configurable)
    DISSIM_WARNING_THRESHOLD <- if (exists("DISSIM_WARNING_THRESHOLD")) {
        DISSIM_WARNING_THRESHOLD 
    } else {
        0.08  # Conservative default
    }
    
    # Calculate degradation status for each study at each level
    study_status <- master_verdicts %>%
        mutate(
            # Use confirmed threshold from consensus voting (v3.2+), or legacy min()
            First_Threshold = pmap_dbl(
                list(
                  Threshold_Observed_Changepoint, Threshold_Observed_Cutoff, Threshold_Observed_NullModel,
                  if ("Confirmed_Threshold_BP" %in% names(.)) Confirmed_Threshold_BP else rep(NA_real_, n()),
                  if ("Consensus_Status"       %in% names(.)) Consensus_Status       else rep(NA_character_, n())
                ),
                compute_first_threshold
            ),

            # Classification logic
            Status = case_when(
                Is_Outlier                                                         ~ "EXCLUDED",
                !is.na(Consensus_Status) & Consensus_Status == "HARD_FAIL"        ~ "FAIL",
                !is.na(Consensus_Status) & Consensus_Status == "WARNING_SINGLE"   ~ "CAUTION",
                is.na(First_Threshold)                                             ~ "PASS",
                First_Threshold >= Threshold_Required                              ~ "PASS",
                First_Threshold >= (Threshold_Required * 0.85)                    ~ "MARGINAL",
                TRUE                                                               ~ "FAIL"
            )
        )
    
    # Check β-diversity at target trim for PASS studies
    # Load all study data to check dissimilarity values
    study_dissim_warnings <- list()
    
    for (study_name in unique(study_status$Study)) {
        # Try to load study data
        study_file_pattern <- paste0("*", study_name, "*results.rds")
        study_files <- list.files(
            path = file.path(OUTDIR, "real_data_results"),
            pattern = study_file_pattern,
            recursive = TRUE,
            full.names = TRUE
        )
        
        if (length(study_files) > 0) {
            tryCatch({
                study_data <- readRDS(study_files[1])
                target_trim <- study_data$target_trim_bp
                
                if (!is.null(study_data$dissim_data) && !is.na(target_trim) && target_trim > 0) {
                    # Check each level
                    for (level in unique(study_status$Level)) {
                        study_level_rows <- study_status %>% 
                            filter(Study == study_name, Level == level, Status == "PASS")
                        
                        if (nrow(study_level_rows) > 0) {
                            dissim_at_target <- study_data$dissim_data %>%
                                filter(Level == level, Trim_BP == target_trim) %>%
                                pull(Dissimilarity)
                            
                            if (length(dissim_at_target) > 0 && !is.na(dissim_at_target[1])) {
                                dissim_value <- dissim_at_target[1]
                                
                                if (dissim_value >= DISSIM_WARNING_THRESHOLD) {
                                    key <- paste(study_name, level, sep = "___")
                                    study_dissim_warnings[[key]] <- dissim_value
                                    
                                    # Update status to PASS_WARNING
                                    study_status <- study_status %>%
                                        mutate(Status = ifelse(
                                            Study == study_name & Level == level & Status == "PASS",
                                            "PASS_WARNING",
                                            Status
                                        ))
                                }
                            }
                        }
                    }
                }
            }, error = function(e) {
                message(paste("  Warning: Could not load study data for", study_name))
            })
        }
    }
    
    # Summarize by taxonomic level
    level_summary <- study_status %>%
        group_by(Level) %>%
        summarise(
            Total_Studies    = n(),
            Passed           = sum(Status %in% c("PASS", "PASS_WARNING"), na.rm = TRUE),
            Passed_Clean     = sum(Status == "PASS",         na.rm = TRUE),
            Passed_Warning   = sum(Status == "PASS_WARNING", na.rm = TRUE),
            Caution_Count    = sum(Status == "CAUTION",      na.rm = TRUE),
            Marginal         = sum(Status == "MARGINAL",     na.rm = TRUE),
            Failed           = sum(Status == "FAIL",         na.rm = TRUE),
            Excluded         = sum(Status == "EXCLUDED",     na.rm = TRUE),

            # Study lists
            Conservative_Studies = paste(Study[Status == "PASS"],    collapse = "; "),
            Conservative_Count   = sum(Status == "PASS"),
            Moderate_Studies     = paste(Study[Status %in% c("PASS", "MARGINAL")], collapse = "; "),
            Moderate_Count       = sum(Status %in% c("PASS", "MARGINAL")),
            Caution_Studies      = paste(Study[Status == "CAUTION"],  collapse = "; "),
            Warning_Studies      = paste(Study[Status == "PASS_WARNING"], collapse = "; "),

            Pass_Rate     = round(Passed_Clean / Total_Studies * 100, 1),
            Moderate_Rate = round((Passed_Clean + Marginal) / Total_Studies * 100, 1),
            Warning_Rate  = round(Passed_Warning / Total_Studies * 100, 1),
            Caution_Rate  = round(Caution_Count  / Total_Studies * 100, 1),
            .groups = "drop"
        ) %>%
        arrange(desc(Moderate_Rate), desc(Pass_Rate))

    # Identify best level
    if (nrow(level_summary) > 0) {
        best_level              <- level_summary$Level[1]
        best_conservative_count <- level_summary$Conservative_Count[1]
        best_moderate_count     <- level_summary$Moderate_Count[1]
        best_caution_count      <- level_summary$Caution_Count[1]
        best_warning_count      <- level_summary$Passed_Warning[1]
        best_pass_rate          <- level_summary$Pass_Rate[1]
        best_moderate_rate      <- level_summary$Moderate_Rate[1]

        use_moderate <- (best_pass_rate < 60 && best_moderate_rate >= 60)

        if (use_moderate) {
            recommended_studies <- level_summary$Moderate_Studies[1]
            recommended_count   <- best_moderate_count
            inclusion_type      <- "Moderate"
        } else {
            recommended_studies <- level_summary$Conservative_Studies[1]
            recommended_count <- best_conservative_count
            inclusion_type <- "Conservative"
        }
        
    } else {
        best_level <- "Family"
        best_pass_rate <- 0
        best_moderate_rate <- 0
        best_warning_count <- 0
        recommended_studies <- "None"
        recommended_count <- 0
        inclusion_type <- "Conservative"
        use_moderate <- FALSE
    }
    
    total_studies <- length(unique(master_verdicts$Study))
    total_analyses <- nrow(master_verdicts)
    
    # Format display table
    display_table <- level_summary %>%
        select(Level, Pass_Rate, Caution_Rate, Warning_Rate, Moderate_Rate,
               Passed_Clean, Caution_Count, Passed_Warning, Marginal, Failed, Excluded) %>%
        rename(
            `Level`          = Level,
            `Pass (%)`       = Pass_Rate,
            `Caution (%)`    = Caution_Rate,
            `Pass+Warn (%)` = Warning_Rate,
            `+Marginal (%)` = Moderate_Rate,
            `n_Pass`         = Passed_Clean,
            `n_Caut`         = Caution_Count,
            `n_Warn`         = Passed_Warning,
            `n_Marg`         = Marginal,
            `n_Fail`         = Failed,
            `n_Excl`         = Excluded
        )
    
    # Build β-diversity warning section
    warning_section <- if (best_warning_count > 0) {
        warning_studies_text <- level_summary %>% 
            filter(Level == best_level) %>% 
            pull(Warning_Studies)
        
        paste(
            "",
            sprintf("⚠ BETA DIVERSITY WARNINGS (%d studies):", best_warning_count),
            sprintf("  Studies with β-div > %.2f at target trim:", DISSIM_WARNING_THRESHOLD),
            paste0("  ", warning_studies_text),
            "",
            "  These studies passed statistical tests but show elevated",
            "  community dissimilarity. Review individual reports carefully.",
            sep = "\n"
        )
    } else {
        ""
    }
    
    # Build recommendation text
    recommendation_text <- paste(
        "═══════════════════════════════════════════════════════════",
        "    MESAP STUDY INCLUSION RECOMMENDATIONS",
        "═══════════════════════════════════════════════════════════",
        "",
        "RECOMMENDED TAXONOMIC LEVEL",
        paste0("  → ", best_level, " (", best_moderate_rate, "% inclusion rate)"),
        "",
        if (use_moderate) {
            paste0("⚠ WARNING: Only ", best_pass_rate, "% of studies fully passed at this level.",
                   "\n  Review marginal studies carefully before inclusion.")
        } else if (best_warning_count > 0) {
            paste0("✓ ", best_pass_rate, "% passed quality checks.",
                   "\n  ⚠ ", best_warning_count, " studies show elevated β-diversity (see below).")
        } else {
            paste0("✓ ", best_pass_rate, "% of studies fully passed quality checks.")
        },
        "",
        "INCLUSION CRITERIA",
        sprintf("  Conservative (n=%d): Studies with NO detected degradation",
                level_summary$Conservative_Count[level_summary$Level == best_level]),
        sprintf("  Moderate (n=%d): Includes marginal cases (within 15%% of target)",
                level_summary$Moderate_Count[level_summary$Level == best_level]),
        sprintf("  Caution (n=%d): Passed but only 1 method triggered (WARNING_SINGLE)",
                level_summary$Caution_Count[level_summary$Level == best_level]),
        sprintf("  Warning (n=%d): Passed but β-diversity > %.2f",
                best_warning_count, DISSIM_WARNING_THRESHOLD),
        "",
        sprintf("RECOMMENDED STUDIES (%s Inclusion):", inclusion_type),
        paste0("  ", recommended_studies),
        warning_section,
        "═══════════════════════════════════════════════════════════",
        "    SUMMARY BY TAXONOMIC LEVEL",
        "═══════════════════════════════════════════════════════════",
        "",
        paste(capture.output(print(display_table, row.names = FALSE)), collapse = "\n"),
        "",
        "═══════════════════════════════════════════════════════════",
        "    ANALYSIS SUMMARY",
        "═══════════════════════════════════════════════════════════",
        paste("  Total Studies Analyzed:", total_studies),
        paste("  Total Level×Study Combinations:", total_analyses),
        sprintf("  β-diversity Warning Threshold: %.2f", DISSIM_WARNING_THRESHOLD),
        paste("  Generated:", Sys.Date()),
        "",
        "═══════════════════════════════════════════════════════════",
        "    CRITICAL: REVIEW INDIVIDUAL REPORTS",
        "═══════════════════════════════════════════════════════════",
        "",
        "⚠ IMPORTANT: Automated threshold detection may be overly strict",
        "  or fail to detect degradation in cases with large initial",
        "  diversity shifts. Individual study reports MUST be reviewed",
        "  before making final inclusion decisions.",
        "",
        "  Key considerations:",
        "  • Verify degradation points align with biological expectations",
        "  • Check for false positives from stochastic variation",
        "  • Assess if 'marginal' studies show real vs. technical effects",
        "  • Review taxon-level impacts for ecological relevance",
        "  • Studies with β-diversity warnings may show gradual drift",
        "",
        "═══════════════════════════════════════════════════════════",
        sep = "\n"
    )
    
    # Create plot
    p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                label = recommendation_text,
                hjust = 0.5, vjust = 0.5, 
                size = 3, family = "mono", lineheight = 1.05) +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = "MESAP Study Inclusion Guide") +
        theme_void() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
            plot.margin = margin(20, 20, 20, 20)
        )
    
    return(p)
}

# ------------------------------------------------------------------------------
# Function: create_method_performance_plot
# Description: Compares which detection methods (Changepoint, Cutoff, Null Model)
#              are triggering the most often and at what thresholds.
# ------------------------------------------------------------------------------
create_method_performance_plot <- function() {
    message("--- Creating method performance comparison plot ---")

    if (nrow(master_verdicts) == 0) {
        message("⚠️ No verdict data available for method comparison")
        return(NULL)
    }

    method_data <- master_verdicts %>%
        dplyr::select(Level, Threshold_Observed_Changepoint, Threshold_Observed_Cutoff, Threshold_Observed_NullModel) %>%
        tidyr::pivot_longer(
            cols = starts_with("Threshold_Observed_"),
            names_to = "Method",
            values_to = "Threshold",
            values_drop_na = TRUE
        ) %>%
        dplyr::mutate(Method = stringr::str_replace(Method, "Threshold_Observed_", "")) %>%
        dplyr::group_by(Level, Method) %>%
        dplyr::summarise(Median_Threshold = median(Threshold, na.rm = TRUE), .groups = "drop") %>%
        tidyr::complete(Level, Method)

    if (nrow(method_data) == 0) {
        message("⚠️ No method data available after processing")
        return(NULL)
    }

    p <- ggplot2::ggplot(method_data, aes(x = Level, y = Median_Threshold, color = Method)) +
        ggplot2::geom_segment(
            aes(xend = Level, yend = 0),
            position = position_dodge(width = 0.7),
            linewidth = 0.7,
            alpha = 0.5
        ) +
        ggplot2::geom_point(
            position = position_dodge(width = 0.7),
            size = 4,
            na.rm = TRUE
        ) +
        ggplot2::geom_text(
            aes(label = round(Median_Threshold)),
            position = position_dodge(width = 0.7),
            vjust = -1.5,
            size = 3,
            na.rm = TRUE,
            fontface = "bold"
        ) +
        ggplot2::scale_color_brewer(palette = "Set2", name = "Detection Method") +
        ggplot2::labs(
            title = "Detection Method Performance Comparison",
            subtitle = "Median degradation threshold detected by each method per taxonomic level",
            x = "Taxonomic Level",
            y = "Median Degradation Threshold (bp)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = "top",
            legend.title = element_text(face = "bold"),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank()
        )

    return(p)
}

# ==============================================================================
# SECTION 6: MAIN PROCESSING LOOP
# This loop iterates through every study's RDS result file, generates the
# individual plots, and compiles them into per-study PDF reports.
# ==============================================================================
check_report_exists <- function(study_name, primer) {
  file.exists(here(output_dir, paste0("Report_Fixed_", study_name, "_", primer, ".pdf")))
}

message("--- Processing studies with file existence checks ---")

real_data_files <- list.files(path = file.path(OUTDIR, "real_data_results"),
                             pattern = "_results\\.rds$",
                             recursive = TRUE,
                             full.names = TRUE)

if (length(real_data_files) == 0) {
    stop("No real data result files found in 'output/real_data_results'.")
}

message(paste("Found", length(real_data_files), "real data files."))

# Log configuration status
if (FORCE_REGENERATE) {
    message("🔄 FORCE_REGENERATE = TRUE: Will regenerate all reports even if they exist.")
} else {
    message("⚡ FORCE_REGENERATE = FALSE: Will skip existing reports to save time.")
}

if (SKIP_INDIVIDUAL) {
    message("⭐️ SKIP_INDIVIDUAL = TRUE: Skipping all individual reports, will only generate master summary.")
}

successful_reports <- 0
skipped_reports <- 0

# ------------------------------------------------------------------------------
# 6A. Individual Report Generation Loop
# ------------------------------------------------------------------------------
if (!SKIP_INDIVIDUAL) {
    for (i in seq_along(real_data_files)) {
        real_file <- real_data_files[i]

        message(paste("\n--- Processing study", i, "of", length(real_data_files), ":", basename(real_file), "---"))

        tryCatch({
            current_study_data <- readRDS(real_file)
            study_name <- current_study_data$study_name

            # --- DEBUGGING FILE VERSION ---
            # Ensures we are not processing stale RDS files from previous pipeline versions
            message(paste("DEBUG: Reading exact file path:", real_file))
            file_info <- file.info(real_file)
            message(paste("DEBUG: File last modified:", file_info$mtime))
            message(paste("DEBUG: File size:", file_info$size, "bytes"))

            if (is.null(current_study_data$taxon_impacts)) {
              message("!!! CRITICAL WARNING: This file lacks 'taxon_impacts'. It is an OLD version.")
              message("!!! ACTION REQUIRED: Re-run Step 06 (analyse_real) for this specific study.")
            } else {
              message("DEBUG: 'taxon_impacts' found. Data is current.")
            }
            # ------------------------------

            primer <- current_study_data$primer_name
            message(paste("  → Loaded study:", study_name, "| Primer:", primer))

            # Check if report already exists
            if (!FORCE_REGENERATE && check_report_exists(study_name, primer)) {
                message("  ⭐️ Report already exists - SKIPPING (set FORCE_REGENERATE=TRUE to regenerate)")
                skipped_reports <- skipped_reports + 1
                rm(current_study_data); gc()
                next
            }

            # -----------------------------------------------------------
            # DETERMINE PLOTTING RANGES (MODE-AWARE)
            # -----------------------------------------------------------

            # 1. Observed Max (BP) - Where the data actually stops
            obs_max_step_bp <- if (!is.null(current_study_data$max_valid_trim_step) && !is.na(current_study_data$max_valid_trim_step)) {
                current_study_data$max_valid_trim_step
            } else {
                if (!is.null(current_study_data$dissim_data) && nrow(current_study_data$dissim_data) > 0) {
                    max(current_study_data$dissim_data$Trim_BP, na.rm = TRUE)
                } else if (!is.null(current_study_data$retention_data) && nrow(current_study_data$retention_data) > 0) {
                    max(current_study_data$retention_data$Trim_BP, na.rm = TRUE)
                } else {
                    300 # Fallback 30 steps * 10bp
                }
            }

            # 2. Axis max: extend to whichever is larger — data endpoint or required
            #    threshold — plus one increment buffer so the green line is always
            #    visible. Never pad to an arbitrary fixed step count, which causes
            #    spurious imputation of data loss on short-amplicon studies.
            actual_increment_study <- if (!is.null(current_study_data$increment) &&
                                           !is.na(current_study_data$increment)) {
                current_study_data$increment
            } else { 100 }

            required_bp_for_study <- tryCatch({
                req <- master_verdicts %>%
                    dplyr::filter(Study == study_name) %>%
                    dplyr::pull(Threshold_Required)
                if (length(req) > 0 && any(!is.na(req))) max(req, na.rm = TRUE) else NA_real_
            }, error = function(e) NA_real_)

            axis_max_bp <- if (!is.na(required_bp_for_study) && required_bp_for_study > obs_max_step_bp) {
                required_bp_for_study + actual_increment_study   # required threshold beyond data — show it
            } else {
                obs_max_step_bp + actual_increment_study          # data covers required — just add one step buffer
            }
            plot_max_step_index <- ceiling(axis_max_bp / actual_increment_study)

            message(paste("  → Axis: data", obs_max_step_bp, "bp | required",
                          ifelse(is.na(required_bp_for_study), "NA", required_bp_for_study),
                          "bp | axis max", plot_max_step_index, "steps"))

            # Filter Verdicts (Mode Aware)
            verdict_data_for_study <- if (ANALYSIS_MODE == "study") {
                master_verdicts %>% filter(Study == study_name)
            } else {
                master_verdicts %>% filter(Primer == primer)
            }

            if (nrow(verdict_data_for_study) == 0) {
                message("     WARNING: No verdict data found for this study. Skipping.")
                rm(current_study_data); gc()
                next
            }

            # Filter Simulation Baselines (Mode Aware)
            baseline_match_id <- if (ANALYSIS_MODE == "study") study_name else primer

            sim_baseline_for_primer <- if (nrow(sim_baselines) > 0) {
                sim_baselines %>% filter(task_id == baseline_match_id)
            } else {
                tibble()
            }

            retention_baseline_for_primer <- if (nrow(sim_retention_curves) > 0) {
                sim_retention_curves %>% filter(task_id == baseline_match_id)
            } else {
                tibble()
            }

            # Create Temp Dir for PNGs
            temp_png_dir <- file.path(output_dir, "temp_pngs", study_name)
            dir.create(temp_png_dir, recursive = TRUE, showWarnings = FALSE)
            png_files_for_pdf <- c()

            # Plot Loop (Phylum -> ASV)
            taxonomic_levels_to_plot <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")
            for (current_level in taxonomic_levels_to_plot) {
                message(paste("  → Generating plots for level:", current_level))

                sim_baseline_level_data <- sim_baseline_for_primer %>% filter(Level == current_level)
                sim_retention_level_data <- retention_baseline_for_primer %>% filter(Level == current_level)
                level_verdict <- verdict_data_for_study %>% filter(Level == current_level)

                if (nrow(level_verdict) == 0) {
                    message(paste("       WARNING: No verdict data for level:", current_level, "- skipping this level."))
                    next
                }

                # --- VERDICT LOGIC (v3.2 consensus-aware) ---

                # 1. Check Outlier Status
                is_outlier <- FALSE
                if ("Is_Outlier" %in% names(level_verdict)) {
                  is_outlier <- any(level_verdict$Is_Outlier, na.rm = TRUE)
                }

                # 2. Extract consensus-aware threshold and method label
                #    get_first_degradation_threshold() now returns NA for WARNING_SINGLE,
                #    so the red line will only appear for CONFIRMED degradation.
                observed_thresh_bp <- get_first_degradation_threshold(level_verdict)
                which_method       <- get_first_degradation_method(level_verdict)
                required_thresh    <- if ("Threshold_Required" %in% names(level_verdict)) level_verdict$Threshold_Required[1] else NA

                # Read consensus status directly for verdict label
                consensus_status <- if ("Consensus_Status" %in% names(level_verdict)) {
                  level_verdict$Consensus_Status[1]
                } else {
                  NA_character_
                }

                # Pre-trim quality flag
                pretrim_flagged <- if ("Pre_Trim_Quality_Flag" %in% names(level_verdict)) {
                  isTRUE(level_verdict$Pre_Trim_Quality_Flag[1])
                } else {
                  FALSE
                }

                # 3. Determine Final Verdict and Display Strings
                if (is_outlier) {
                  mainverdict        <- "FAIL (Outlier)"
                  subtitle_obs       <- "NA (Alignment Mismatch)"
                  observed_thresh_bp <- NA
                  warning_thresh_bp  <- NA_real_
                  required_thresh    <- NA
                  level_verdict_for_plotting <- level_verdict %>%
                    mutate(
                      Threshold_Required             = NA_real_,
                      Threshold_Observed_Changepoint = NA_real_,
                      Threshold_Observed_Cutoff      = NA_real_,
                      Threshold_Observed_NullModel   = NA_real_,
                      Confirmed_Threshold_BP         = NA_real_
                    )

                } else {
                  level_verdict_for_plotting <- level_verdict

                  mainverdict <- if (!is.na(consensus_status)) {
                    # v3.2 path: read status directly
                    if (consensus_status == "HARD_FAIL") {
                      "FAIL (quality gate)"   # BC ceiling or retention floor breached
                    } else if (consensus_status == "CONFIRMED") {
                      if (!is.na(observed_thresh_bp) && !is.na(required_thresh) &&
                          observed_thresh_bp < required_thresh) {
                        "FAIL"
                      } else {
                        "PASS"
                      }
                    } else if (consensus_status == "WARNING_SINGLE") {
                      "CAUTION"   # single-method trigger — not confirmed degradation
                    } else {
                      "PASS"      # NONE
                    }
                  } else {
                    # Legacy fallback
                    if (is.na(observed_thresh_bp) || is.na(required_thresh)) {
                      "PASS"
                    } else if (observed_thresh_bp < required_thresh) {
                      "FAIL"
                    } else {
                      "PASS"
                    }
                  }

                  # Compute warning_thresh_bp for WARNING_SINGLE — used for orange line on plots
                  warning_thresh_bp <- NA_real_
                  if (!is.na(consensus_status) && consensus_status == "WARNING_SINGLE") {
                    raw_min <- min(c(level_verdict$Threshold_Observed_Changepoint,
                                     level_verdict$Threshold_Observed_Cutoff,
                                     level_verdict$Threshold_Observed_NullModel), na.rm = TRUE)
                    if (is.finite(raw_min)) warning_thresh_bp <- raw_min
                  }

                  subtitle_obs <- if (is.na(observed_thresh_bp)) {
                    if (!is.na(warning_thresh_bp)) {
                      paste0(warning_thresh_bp, " bp (1-method warning)")
                    } else {
                      "None"
                    }
                  } else {
                    paste(observed_thresh_bp, "bp")
                  }

                  # Append pre-trim flag if set
                  if (pretrim_flagged) subtitle_obs <- paste0(subtitle_obs, " ⚠PreTrim")
                }

                subtitle_txt <- paste("Level:", current_level,
                                      "| Method:", which_method,
                                      "| Observed:", subtitle_obs,
                                      "| Verdict:", mainverdict)

                # Generate the 4 Panels
	p1 <- plot_dissimilarity_robust(current_study_data, sim_baseline_level_data, level_verdict_for_plotting, obs_max_step_bp, current_level, observed_thresh_bp, plot_max_step_index, warning_thresh_bp = warning_thresh_bp)
	p2 <- plot_retention_robust(current_study_data, sim_retention_level_data, level_verdict_for_plotting, obs_max_step_bp, current_level, observed_thresh_bp, plot_max_step_index, warning_thresh_bp = warning_thresh_bp)
	p3 <- plot_taxon_impact_combined(current_study_data, level_verdict_for_plotting, current_level, obs_max_step_bp, observed_thresh_bp, plot_max_step_index, warning_thresh_bp = warning_thresh_bp)
	p4 <- plot_core_taxa_robust(current_study_data, level_verdict_for_plotting, obs_max_step_bp, current_level, observed_thresh_bp, plot_max_step_index, warning_thresh_bp = warning_thresh_bp)

                # Combine into Grid and Save PNG
                tryCatch({
                    plot_grid <- wrap_plots(p1, p2, p3, p4, nrow = 2)
                    annotated_plot_grid <- plot_grid +
                        plot_annotation(
                            title = paste("MESAP Analysis Report:", study_name),
                            subtitle = subtitle_txt,
                            caption = paste("Generated:", Sys.Date(), "| Primer:", primer, "| Valid data to", obs_max_step_bp, "bp"),
                            theme = theme(plot.title = element_text(size = 18, face = "bold"))
                        )

                    temp_png_filename <- file.path(temp_png_dir, paste0("page_", match(current_level, taxonomic_levels_to_plot), ".png"))
                    ggsave(temp_png_filename, annotated_plot_grid, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = "white")
                    png_files_for_pdf <- c(png_files_for_pdf, temp_png_filename)
                    message(paste("    ✅ Created page for level:", current_level))
                    rm(p1, p2, p3, p4, plot_grid, annotated_plot_grid)
                }, error = function(e) {
                    message(paste("    ⛔ Error creating plot grid for level", current_level, ":", conditionMessage(e)))
                })
                gc()
            }

            # Stitch PNGs into PDF using Magick
            if (length(png_files_for_pdf) > 0) {
                message(paste("  → Creating PDF with", length(png_files_for_pdf), "pages..."))
                tryCatch({
                    image_stack <- image_read(png_files_for_pdf)
                    pdf_filename <- here(output_dir, paste0("Report_Fixed_", study_name, "_", primer, ".pdf"))
                    image_write(image_stack, path = pdf_filename, format = "pdf")
                    message(paste("  ✅ Created:", pdf_filename))
                    successful_reports <- successful_reports + 1
                    rm(image_stack)
                }, error = function(e) {
                    message(paste("  ⛔ Error creating PDF:", conditionMessage(e)))
                })
            }

            # Cleanup
            unlink(temp_png_dir, recursive = TRUE)
            rm(current_study_data, verdict_data_for_study, sim_baseline_for_primer, retention_baseline_for_primer)
            gc()

        }, error = function(e) {
            message(paste("  ⛔ FATAL ERROR processing study:", basename(real_file), "->", conditionMessage(e)))
            if (exists("current_study_data")) rm(current_study_data)
            gc()
        })
    }
} else {
    message("⭐️ Skipping all individual report generation because SKIP_INDIVIDUAL is TRUE.")
}

# ------------------------------------------------------------------------------
# 6B. Master Summary Report Generation
# ------------------------------------------------------------------------------
message("\n--- Creating Improved Master Summary Report ---")

tryCatch({
    taxonomic_levels_to_plot <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")
    master_plots_list <- list()

    # Generate Alignment Plots (Map)
    for (level in taxonomic_levels_to_plot) {
        message(paste("  -> Generating alignment plot for level:", level))
        alignment_plot <- tryCatch({
            create_amplicon_alignment_plot(master_manifest = master_manifest, master_verdicts = master_verdicts, taxonomic_level = level)
        }, error = function(e) {
            message(paste("    [WARN] Error creating alignment plot for", level, ":", conditionMessage(e)))
            NULL
        })

        if (!is.null(alignment_plot)) {
            master_plots_list[[paste0("alignment_", level)]] <- alignment_plot
            message(paste("    [OK] Created alignment plot for", level))
        } else {
            message(paste("    [SKIP] Skipping alignment plot for", level, "(no data or error)"))
        }
    }

    # Generate Forest Plot
    message("  -> Generating trimming recommendations...")
    master_plots_list$recommendations <- tryCatch({
        create_trimming_recommendations()
    }, error = function(e) {
        message(paste("    [WARN] Error creating recommendations:", conditionMessage(e)))
        NULL
    })

    # Generate Cheat Sheet
    message("  -> Generating practical guide...")
    master_plots_list$practical_guide <- tryCatch({
        create_practical_guide()
    }, error = function(e) {
        message(paste("    [WARN] Error creating practical guide:", conditionMessage(e)))
        NULL
    })

    # Generate Method Comparison
    message("  -> Generating method performance plot...")
    master_plots_list$method_performance <- tryCatch({
        create_method_performance_plot()
    }, error = function(e) {
        message(paste("    [WARN] Error creating method performance plot:", conditionMessage(e)))
        NULL
    })

    # Remove NULLs
    master_plots_list <- purrr::compact(master_plots_list)

    # Stitch Master PDF
    if (length(master_plots_list) > 0) {
        message(paste("  -> Creating PDF with", length(master_plots_list), "pages..."))

        temp_master_dir <- file.path(output_dir, "temp_master")
        dir.create(temp_master_dir, recursive = TRUE, showWarnings = FALSE)
        master_png_files <- c()

        for (i in seq_along(master_plots_list)) {
            plot_name <- names(master_plots_list)[i]
            current_plot <- master_plots_list[[i]]

            png_filename <- file.path(temp_master_dir, paste0(sprintf("%02d", i), "_", plot_name, ".png"))

            tryCatch({
                ggsave(png_filename, current_plot,
                       width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = "white")
                master_png_files <- c(master_png_files, png_filename)
                message(paste("    [OK] Saved page", i, ":", plot_name))
            }, error = function(e) {
                message(paste("    [ERR] Error saving", plot_name, ":", conditionMessage(e)))
            })
        }

        if (length(master_png_files) > 0) {
            message("  -> Combining pages into Master Summary PDF...")
            tryCatch({
                master_images <- image_read(sort(master_png_files))
                master_pdf_filename <- here(output_dir, "MESAP_Master_Summary_Report.pdf")
                image_write(master_images, path = master_pdf_filename, format = "pdf")
                message(paste("  [OK] Created Master Summary:", master_pdf_filename))
            }, error = function(e) {
                message(paste("  [ERR] Error creating master PDF:", conditionMessage(e)))
            })
        } else {
            message("  [WARN] No master summary pages were successfully created.")
        }

        unlink(temp_master_dir, recursive = TRUE)

    } else {
        message("  [WARN] No master summary plots could be created - all plot generation failed.")
    }

}, error = function(e) {
    message(paste("  [FATAL] FATAL ERROR creating improved master summary:", conditionMessage(e)))
    message(paste("     Traceback:", paste(capture.output(traceback()), collapse = "\n")))
})

# Final Stats Log
message("\n=== REPORT GENERATION COMPLETE ===")
if(exists("successful_reports") && exists("skipped_reports")) {
    message(paste("Successful individual reports:", successful_reports))
    message(paste("Skipped individual reports:", skipped_reports))
    message(paste("Total processed:", successful_reports + skipped_reports))
}
message(paste("Configuration used: PLOT_WIDTH =", PLOT_WIDTH, "| PLOT_HEIGHT =", PLOT_HEIGHT, "| DPI =", PLOT_DPI))
