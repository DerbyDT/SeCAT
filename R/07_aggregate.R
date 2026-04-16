#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT:   scripts/07_aggregate.R
# PIPELINE: SeCAT (Sequence Consensus Analysis Tool)
# PHASE:    Phase 5: Aggregation & Statistical Verdicts
# VERSION:  3.0 (Memory Optimized & Dual-Mode)
# AUTHOR:   [Author Name]
#
# PURPOSE:
#   This is the "Reducer" script. It collects the results from:
#   1. All Real Data Analysis jobs (Phase 4)
#   2. All Simulation jobs (Phase 3)
#
#   It performs the final statistical comparisons to generate "Verdicts" for
#   each study:
#   - Is the observed degradation significantly worse than the null model?
#   - Where is the "tipping point" (changepoint)?
#   - Does the study fail the absolute cutoff?
#
#   Finally, it produces the Master CSV tables used by the Report Generator.
#
# MEMORY OPTIMIZATION:
#   Because simulation data can be massive (TB scale), this script uses
#   batch processing and streaming writes to keep RAM usage low (<10GB).
# ==============================================================================

log_and_flush <- function(message) {
  cat(paste(Sys.time(), "|", message, "\n"))
  flush.console()
}

# ==============================================================================
# SECTION 1: STATISTICAL HELPER FUNCTIONS
# ==============================================================================

# Function: get_first_degradation_threshold
# Description: Helper to find the most conservative (lowest BP) threshold triggered.
get_first_degradation_threshold <- function(vdf) {
  all_thresholds <- c(vdf$Threshold_Observed_Changepoint, vdf$Threshold_Observed_Cutoff, vdf$Threshold_Observed_NullModel)
  all_thresholds <- all_thresholds[!is.na(all_thresholds)]
  if (length(all_thresholds) == 0) return(NA_real_) else return(min(all_thresholds))
}

get_levels <- function() {
  return(c("Phylum", "Class", "Order", "Family", "Genus", "ASV"))
}

# ------------------------------------------------------------------------------
# Function: select_changepoint_penalty_cv
# Description: Selects the optimal penalty for PELT changepoint detection using
#              bootstrap cross-validation. This prevents "over-fitting" where
#              tiny fluctuations are flagged as changepoints.
#
# Args:
#   y: Numeric vector of dissimilarity values.
#   trim_bp: Numeric vector of trim positions.
#   penalty_method: "AIC", "SIC", or "MANUAL".
#   scan: Vector of penalty multipliers to test.
#   nboot: Number of bootstrap replicates.
# ------------------------------------------------------------------------------
select_changepoint_penalty_cv <- function(y, trim_bp, penalty_method, scan, nboot, min_trim_bp = 10) {
  requireNamespace("changepoint")

  # Filter valid indices first to avoid repeated subsetting
  valid_indices <- trim_bp > min_trim_bp
  y_filtered <- y[valid_indices]
  trim_bp_filtered <- trim_bp[valid_indices]

  n <- length(y_filtered)
  # Early exit if not enough data
  if(n < 3) return(list(penalty_mult=NA, changepoint_bp=NA, consensus_frac=NA))

  results <- list()

  # Define penalty function once based on method
  penalty_func <- switch(penalty_method,
                         "AIC" = function(mult) mult * 2,
                         "SIC" = function(mult) mult * log(n),
                         "MANUAL" = function(mult) mult * log(n),
                         function(mult) mult * log(n))

  # Pre-calculate penalty values to avoid function calls in loop
  penalty_values <- sapply(scan, penalty_func)

  for (i in seq_along(scan)) {
    mult <- scan[i]
    p_val <- penalty_values[i]

    cp_results <- numeric(nboot)

    # Run bootstrap
    for (b in seq_len(nboot)) {
      # Sample indices with replacement
      idx <- sample.int(n, size = n, replace = TRUE)
      y_boot <- y_filtered[idx]
      tb_boot <- trim_bp_filtered[idx]

      cp_loc <- tryCatch({
        # changepoint::cpt.meanvar is lightweight, but we wrap in tryCatch
        cp <- changepoint::cpt.meanvar(y_boot, method="PELT", penalty="Manual", pen.value=p_val)
        cpts <- changepoint::cpts(cp)
        # Return the BP location of the first changepoint
        if (length(cpts) > 0 && cpts[1] < length(tb_boot)) tb_boot[cpts[1]] else NA_real_
      }, error=function(e) NA_real_)

      cp_results[b] <- cp_loc
    }

    # Calculate consensus statistics
    valid_cps <- cp_results[!is.na(cp_results)]
    if (length(valid_cps) == 0) {
      cp_mode <- NA_real_
      frac_mode <- NA_real_
    } else {
      # Efficient mode calculation
      tab <- table(valid_cps)
      cp_mode <- as.numeric(names(tab)[which.max(tab)])
      # Fraction of runs that agreed on this changepoint (Stability Score)
      frac_mode <- mean(cp_results == cp_mode, na.rm = TRUE)
    }

    results[[i]] <- list(penalty_mult = mult, consensus_cp = cp_mode, consensus_frac = frac_mode)
  }

  penalty_tbl <- dplyr::bind_rows(results) %>% dplyr::filter(!is.na(consensus_cp))

  if (nrow(penalty_tbl) == 0) return(list(penalty_mult=NA, changepoint_bp=NA, consensus_frac=NA))

  # Select the penalty that yielded the most STABLE changepoint (highest consensus fraction)
  max_consensus_frac <- max(penalty_tbl$consensus_frac, na.rm = TRUE)
  best <- penalty_tbl %>% dplyr::filter(consensus_frac == max_consensus_frac) %>% dplyr::arrange(penalty_mult) %>% dplyr::slice(1)

  list(penalty_mult = best$penalty_mult, changepoint_bp = best$consensus_cp, consensus_frac = best$consensus_frac)
}

# ------------------------------------------------------------------------------
# Function: find_degradation_point (Absolute Threshold)
# Description: Identifies where dissimilarity exceeds a fixed cutoff (e.g., 0.15).
#              Also calculates an empirical p-value against the null model at that point.
# ------------------------------------------------------------------------------
find_degradation_point <- function(real_curve, sim_dissim_data, sim_dissim_baseline, level, task_id, distance_cutoff = 0.15, min_trim_bp = 10) {
  sub_real <- real_curve %>% dplyr::filter(Level == level, Trim_BP >= 0) %>% dplyr::arrange(Trim_BP)
  if (nrow(sub_real) < 2) return(list(degradation_bp = NA_real_, empirical_p_value = NA_real_))

  # Calculate rate of change (delta)
  sub_real <- sub_real %>% mutate(delta_D = Dissimilarity - lag(Dissimilarity, default = Dissimilarity[1]))
  steps <- sub_real %>% dplyr::filter(Trim_BP > min_trim_bp)

  for (i in seq_len(nrow(steps))) {
    trim_bp  <- steps$Trim_BP[i]; distance_real <- steps$delta_D[i]
    prev_trim_bp <- sub_real %>% filter(Trim_BP < trim_bp) %>% pull(Trim_BP) %>% max()

    # Get simulation data for this exact step to calculate p-value
    sim_data_for_step <- sim_dissim_data %>% dplyr::filter(task_id == !!task_id, Level == level, Trim_BP %in% c(prev_trim_bp, trim_bp))
    if (nrow(sim_data_for_step) < 2) next

    sim_deltas <- sim_data_for_step %>% group_by(simulation_id) %>% arrange(Trim_BP) %>% summarise(delta_D = Dissimilarity[2] - Dissimilarity[1], .groups = "drop") %>% filter(!is.na(delta_D)) %>% pull(delta_D)

    # Empirical P-Value: (Number of sims >= real + 1) / (Total sims + 1)
    empirical_p <- if(length(sim_deltas) > 0) { (sum(sim_deltas >= distance_real, na.rm = TRUE) + 1) / (length(sim_deltas) + 1) } else { NA_real_ }

    # Trigger condition: Dissimilarity exceeds cutoff
    if (!is.na(distance_real) && distance_real > distance_cutoff) { return(list(degradation_bp = trim_bp, empirical_p_value = empirical_p)) }
  }
  list(degradation_bp = NA_real_, empirical_p_value = NA_real_)
}

# ------------------------------------------------------------------------------
# Function: find_null_model_point (Statistical Threshold)
# Description: Identifies where the real curve deviates significantly (p < 0.05)
#              from the null model distribution for 3 consecutive steps.
# ------------------------------------------------------------------------------
find_null_model_point <- function(real_curve, sim_dissim_data, level, task_id, empirical_p_threshold = 0.05, min_consecutive_steps = 3, min_trim_bp = 10) {
    sub_real <- real_curve %>% dplyr::filter(Level == level, Trim_BP >= 0) %>% dplyr::arrange(Trim_BP)
    if (nrow(sub_real) < 2) return(list(nullmodel_bp = NA_real_, empirical_p_value = NA_real_))

    sub_real <- sub_real %>% mutate(delta_D = Dissimilarity - lag(Dissimilarity, default = Dissimilarity[1]))
    steps <- sub_real %>% dplyr::filter(Trim_BP > min_trim_bp)
    if (nrow(steps) == 0) return(list(nullmodel_bp = NA_real_, empirical_p_value = NA_real_))

    p_values <- numeric(nrow(steps)); names(p_values) <- steps$Trim_BP

    for (i in seq_len(nrow(steps))) {
        trim_bp  <- steps$Trim_BP[i]; real_delta_D <- steps$delta_D[i]
        prev_trim_bp <- sub_real %>% filter(Trim_BP < trim_bp) %>% pull(Trim_BP) %>% max()

        sim_data_for_step <- sim_dissim_data %>% dplyr::filter(task_id == !!task_id, Level == level, Trim_BP %in% c(prev_trim_bp, trim_bp))
        if (nrow(sim_data_for_step) == 0) { p_values[i] <- NA_real_; next }

        sim_deltas <- sim_data_for_step %>% group_by(simulation_id) %>% arrange(Trim_BP) %>% summarise(delta_D = Dissimilarity[2] - Dissimilarity[1], .groups = "drop") %>% filter(!is.na(delta_D)) %>% pull(delta_D)

        p_values[i] <- if(length(sim_deltas) > 0) { (sum(sim_deltas >= real_delta_D, na.rm = TRUE) + 1) / (length(sim_deltas) + 1) } else { NA_real_ }
    }

    # Requirement: 3 consecutive significant steps to avoid noise triggers
    if (length(p_values) >= min_consecutive_steps) {
        for (i in 1:(length(p_values) - min_consecutive_steps + 1)) {
            window <- p_values[i:(i + min_consecutive_steps - 1)]
            if (all(!is.na(window)) && all(window < empirical_p_threshold)) {
                return(list(nullmodel_bp = as.numeric(names(p_values)[i]), empirical_p_value = p_values[i]))
            }
        }
    }
    list(nullmodel_bp = NA_real_, empirical_p_value = NA_real_)
}

# ------------------------------------------------------------------------------
# Function: find_bc_ceiling_point
# Description: Flags when the absolute Bray-Curtis dissimilarity of the REAL
#              data curve first exceeds bc_ceiling at any trim step up to
#              target_trim.  No simulation comparison needed — this is a direct
#              exceedance check.  Returns the first trim step that exceeds the
#              ceiling, or NA if the ceiling is never reached.
# ------------------------------------------------------------------------------
find_bc_ceiling_point <- function(real_curve, level, target_trim,
                                   bc_ceiling = 0.20, min_trim_bp = 50) {
  sub <- real_curve %>%
    dplyr::filter(Level == level, Trim_BP >= min_trim_bp,
                  !is.na(Dissimilarity)) %>%
    dplyr::arrange(Trim_BP)

  if (nrow(sub) == 0) return(list(bc_ceiling_bp = NA_real_))

  # Only test steps up to (and including) target_trim
  if (!is.na(target_trim)) sub <- sub %>% dplyr::filter(Trim_BP <= target_trim)
  if (nrow(sub) == 0) return(list(bc_ceiling_bp = NA_real_))

  hit <- sub %>% dplyr::filter(Dissimilarity > bc_ceiling)
  if (nrow(hit) == 0) return(list(bc_ceiling_bp = NA_real_))

  list(bc_ceiling_bp = hit$Trim_BP[1])
}

# ------------------------------------------------------------------------------
# Function: find_retention_floor_point
# Description: Flags when taxon RETENTION drops below retention_floor at any
#              trim step up to target_trim.  Operates on the retention_data
#              component of the study RDS.  Returns the first trim step that
#              breaches the floor, or NA if retention never falls that low.
# ------------------------------------------------------------------------------
find_retention_floor_point <- function(retention_data, level, target_trim,
                                        retention_floor = 0.70, min_trim_bp = 50) {
  if (is.null(retention_data) || nrow(retention_data) == 0)
    return(list(retention_floor_bp = NA_real_))

  sub <- retention_data %>%
    dplyr::filter(Level == level, Trim_BP >= min_trim_bp,
                  !is.na(Retention)) %>%
    dplyr::arrange(Trim_BP)

  if (nrow(sub) == 0) return(list(retention_floor_bp = NA_real_))

  if (!is.na(target_trim)) sub <- sub %>% dplyr::filter(Trim_BP <= target_trim)
  if (nrow(sub) == 0) return(list(retention_floor_bp = NA_real_))

  # Retention is stored as a proportion (0-1) or percentage (0-100) — normalise
  if (max(sub$Retention, na.rm = TRUE) > 1.5) sub <- sub %>% dplyr::mutate(Retention = Retention / 100)

  hit <- sub %>% dplyr::filter(Retention < retention_floor)
  if (nrow(hit) == 0) return(list(retention_floor_bp = NA_real_))

  list(retention_floor_bp = hit$Trim_BP[1])
}

# ------------------------------------------------------------------------------
# Function: get_consensus_threshold
# Description: Applies consensus voting across ALL detection methods.
#              A degradation is CONFIRMED only when >= 2 methods agree within
#              tolerance_bp of each other.  A single trigger is WARNING_SINGLE.
#              Returns a list with: status, confirmed_bp, n_triggered, n_agree.
# ------------------------------------------------------------------------------
get_consensus_threshold <- function(cp_bp, cut_bp, null_bp,
                                     tolerance_bp = 500) {
  # Collect non-NA thresholds from the three statistical detection methods.
  # BC ceiling and retention floor are hard cutoffs handled separately —
  # they do not participate in consensus voting.
  all_thresh <- c(
    Changepoint = cp_bp,
    Cutoff      = cut_bp,
    NullModel   = null_bp
  )
  triggered <- all_thresh[!is.na(all_thresh)]
  n_trig    <- length(triggered)

  if (n_trig == 0) return(list(status       = "NONE",
                                confirmed_bp  = NA_real_,
                                outer_bp      = NA_real_,
                                inner_bp      = NA_real_,
                                n_triggered   = 0L,
                                n_agree       = 0L))

  if (n_trig == 1) return(list(status       = "WARNING_SINGLE",
                                confirmed_bp  = NA_real_,
                                outer_bp      = triggered[[1]],   # only trigger = both outer and inner
                                inner_bp      = triggered[[1]],
                                n_triggered   = 1L,
                                n_agree       = 1L))

  # Find the pair with the most agreements within tolerance_bp
  best_anchor  <- NA_real_
  best_n_agree <- 0L

  for (anchor in triggered) {
    n_agree <- sum(abs(triggered - anchor) <= tolerance_bp)
    if (n_agree > best_n_agree) {
      best_n_agree <- n_agree
      best_anchor  <- anchor
    }
  }

  if (best_n_agree >= 2) {
    agreeing     <- triggered[abs(triggered - best_anchor) <= tolerance_bp]
    confirmed_bp <- median(agreeing)
    # outer_bp = most lenient (highest) trigger across ALL fired methods
    # inner_bp = most conservative (lowest) trigger across ALL fired methods
    # These define the two-stage caution → red visual bands on the alignment map.
    return(list(status       = "CONFIRMED",
                confirmed_bp  = confirmed_bp,
                outer_bp      = max(triggered),
                inner_bp      = min(triggered),
                n_triggered   = n_trig,
                n_agree       = best_n_agree))
  }

  # All triggers fired but none within tolerance — still record outer/inner for plotting
  return(list(status       = "WARNING_SINGLE",
              confirmed_bp  = NA_real_,
              outer_bp      = max(triggered),
              inner_bp      = min(triggered),
              n_triggered   = n_trig,
              n_agree       = 1L))
}

# ==============================================================================
# SECTION 2: MEMORY-EFFICIENT PROCESSING
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: process_task_simulations_memory_safe
# Description: Loads simulation replicates for a task, interpolates them to the
#              Real Data grid (handling scaled vs absolute mismatches), and
#              computes summary stats. Uses batch loading to save RAM.
# ------------------------------------------------------------------------------
process_task_simulations_memory_safe <- function(current_task_id,
                                                 simulation_results_dir,
                                                 output_dir,
                                                 real_data_results_dir) {
  log_and_flush(paste("  -> Processing simulations for task:", current_task_id))
  task_sim_dir <- file.path(simulation_results_dir, current_task_id)

  if (!dir.exists(task_sim_dir)) {
    log_and_flush(paste("     WARNING: No simulation directory found for", current_task_id))
    return(NULL)
  }

  # ===== LOAD REAL DATA TO GET INCREMENT AND TRIM_BP GRID =====
  real_result_path <- file.path(real_data_results_dir,
                                current_task_id,
                                paste0(current_task_id, "_results.rds"))
  real_increment <- 10  # Default fallback
  real_trim_bp_grid <- NULL

  if (file.exists(real_result_path)) {
    tryCatch({
      real_data <- readRDS(real_result_path)

      # Extract increment
      if (!is.null(real_data$increment) && !is.na(real_data$increment)) {
        real_increment <- real_data$increment
      }

      # Extract actual Trim_BP grid from dissim_data
      if (!is.null(real_data$dissim_data) && nrow(real_data$dissim_data) > 0) {
        real_trim_bp_grid <- sort(unique(real_data$dissim_data$Trim_BP))
        log_and_flush(paste("     Real data uses increment =", real_increment, "bp"))
        log_and_flush(paste(
          "     Real Trim_BP grid:",
          paste(head(real_trim_bp_grid, 5), collapse = ", "),
          "...",
          max(real_trim_bp_grid)
        ))
      } else {
        log_and_flush("     WARNING: No dissim_data in real data. Using fallback grid.")
        # Fallback: generate grid from 0 to max_valid_trim_step
        max_bp <- if (!is.null(real_data$max_valid_trim_step)) real_data$max_valid_trim_step else 11500
        real_trim_bp_grid <- seq(0, max_bp, by = real_increment)
      }

      rm(real_data)
      gc()

    }, error = function(e) {
      log_and_flush(paste(
        "     WARNING: Could not load real data for",
        current_task_id, "- Error:", conditionMessage(e)
      ))
      log_and_flush(paste("     Using default increment =", real_increment, "bp"))
    })
  } else {
    log_and_flush(paste("     WARNING: No real data file found at", real_result_path))
    log_and_flush(paste("     Using default increment =", real_increment, "bp"))
  }

  # If we still don't have a grid, create a sensible default
  if (is.null(real_trim_bp_grid)) {
    real_trim_bp_grid <- seq(0, 11500, by = real_increment)
    log_and_flush(paste("     Using default Trim_BP grid: 0 to 11500 by", real_increment))
  }

  # ===== LOAD AND BATCH PROCESS SIMULATION FILES =====
  seed_dirs <- list.dirs(task_sim_dir, recursive = FALSE, full.names = TRUE)
  seed_files <- file.path(seed_dirs, "results.rds")
  seed_files <- seed_files[file.exists(seed_files)]

  log_and_flush(paste("     Found", length(seed_files), "simulation files for", current_task_id))
  if (length(seed_files) == 0) return(NULL)

  batch_size <- 25
  all_dissim_data <- list()
  all_retention_data <- list()

  for (batch_start in seq(1, length(seed_files), by = batch_size)) {
    batch_end   <- min(batch_start + batch_size - 1, length(seed_files))
    batch_files <- seed_files[batch_start:batch_end]
    log_and_flush(paste(
      "     Processing batch",
      ceiling(batch_start / batch_size),
      "of",
      ceiling(length(seed_files) / batch_size)
    ))

    batch_dissim    <- list()
    batch_retention <- list()

    for (sim_file in batch_files) {

      aligned_dissim    <- NULL  # ensure defined each iteration
      aligned_retention <- NULL

      tryCatch({
        sim_result <- readRDS(sim_file)
        if (is.null(sim_result$task_id)) sim_result$task_id <- current_task_id
        simulation_id <- paste(sim_result$task_id, sim_result$seed, sep = "_")

        # ===== INTERPOLATE DISSIMILARITY DATA TO REAL GRID =====
        if (!is.null(sim_result$dissim_data) && nrow(sim_result$dissim_data) > 0) {

          aligned_dissim <- sim_result$dissim_data %>%
            dplyr::group_by(Level) %>%
            dplyr::arrange(Trim_BP) %>%
            dplyr::summarise(
              aligned_data = list({
                # Capture original x/y vectors BEFORE defining new columns
                x_orig <- Trim_BP
                y_orig <- Dissimilarity

                tibble::tibble(
                  Trim_BP = real_trim_bp_grid,
                  Dissimilarity = approx(
                    x = x_orig,
                    y = y_orig,
                    xout = real_trim_bp_grid,
                    rule = 2
                  )$y
                )
              }),
              .groups = "drop"
            ) %>%
            tidyr::unnest(aligned_data) %>%
            dplyr::mutate(
              task_id      = sim_result$task_id,
              seed         = sim_result$seed,
              simulation_id = simulation_id,
              mode         = "both",
              TrimStep     = Trim_BP / real_increment
            )
        }

        # ===== INTERPOLATE RETENTION DATA TO REAL GRID =====
        if (!is.null(sim_result$retention_data) && nrow(sim_result$retention_data) > 0) {

          aligned_retention <- sim_result$retention_data %>%
            dplyr::group_by(Level) %>%
            dplyr::arrange(Trim_BP) %>%
            dplyr::summarise(
              aligned_data = list({
                # Capture original x/y vectors BEFORE defining new columns
                x_orig <- Trim_BP
                y_orig <- Retention

                tibble::tibble(
                  Trim_BP = real_trim_bp_grid,
                  Retention = approx(
                    x = x_orig,
                    y = y_orig,
                    xout = real_trim_bp_grid,
                    rule = 2
                  )$y
                )
              }),
              .groups = "drop"
            ) %>%
            tidyr::unnest(aligned_data) %>%
            dplyr::mutate(
              task_id      = sim_result$task_id,
              seed         = sim_result$seed,
              simulation_id = simulation_id,
              mode         = "both",
              TrimStep     = Trim_BP / real_increment
            )
        }

        # Append only if successfully created
        if (!is.null(aligned_dissim)) {
          batch_dissim[[length(batch_dissim) + 1]] <- aligned_dissim
        }
        if (!is.null(aligned_retention)) {
          batch_retention[[length(batch_retention) + 1]] <- aligned_retention
        }

      }, error = function(e) {
        log_and_flush(paste(
          "       WARNING: Failed to process",
          sim_file, "- Error:", conditionMessage(e)
        ))
      })
    }

    if (length(batch_dissim) > 0) {
      all_dissim_data[[length(all_dissim_data) + 1]] <- dplyr::bind_rows(batch_dissim)
    }
    if (length(batch_retention) > 0) {
      all_retention_data[[length(all_retention_data) + 1]] <- dplyr::bind_rows(batch_retention)
    }
    rm(batch_dissim, batch_retention)
    gc()
  }

  # ===== CALCULATE BASELINE STATISTICS ON ALIGNED DATA =====
  task_results <- list()
  if (length(all_dissim_data) > 0) {
    log_and_flush(paste("     Calculating baseline statistics for", current_task_id))
    task_dissim <- dplyr::bind_rows(all_dissim_data)

    task_baseline <- task_dissim %>%
      dplyr::group_by(task_id, Level, Trim_BP, TrimStep) %>%
      dplyr::summarise(
        n_simulations     = dplyr::n(),
        MeanDissimilarity = mean(Dissimilarity, na.rm = TRUE),
        SD_Dissimilarity  = sd(Dissimilarity, na.rm = TRUE),
        CILower           = stats::quantile(Dissimilarity, 0.025, na.rm = TRUE),
        CIUpper           = stats::quantile(Dissimilarity, 0.975, na.rm = TRUE),
        .groups           = "drop"
      )

    task_baseline_file <- file.path(output_dir,
                                    paste0("simulation_baseline_", current_task_id, ".csv"))
    readr::write_csv(task_baseline, task_baseline_file)
    log_and_flush(paste("     ✅ Created", task_baseline_file))

    task_results$baseline_stats   <- task_baseline
    task_results$full_dissim_data <- task_dissim

    if (length(all_retention_data) > 0) {
      task_retention <- dplyr::bind_rows(all_retention_data)

      task_retention_summary <- task_retention %>%
        dplyr::group_by(task_id, Level, Trim_BP, TrimStep) %>%
        dplyr::summarise(
          MeanRetention = mean(Retention, na.rm = TRUE),
          SD            = sd(Retention, na.rm = TRUE),
          CILower       = stats::quantile(Retention, 0.025, na.rm = TRUE),
          CIUpper       = stats::quantile(Retention, 0.975, na.rm = TRUE),
          .groups       = "drop"
        )

      task_retention_file <- file.path(output_dir,
                                       paste0("simulation_retention_", current_task_id, ".csv"))
      readr::write_csv(task_retention_summary, task_retention_file)
      log_and_flush(paste("     ✅ Created", task_retention_file))

      task_results$retention_data <- task_retention
    }
  }

  rm(all_dissim_data, all_retention_data)
  if (exists("task_dissim"))   rm(task_dissim)
  if (exists("task_retention")) rm(task_retention)
  gc()

  return(task_results)
}

# ==============================================================================
# SECTION 3: PER-STUDY TABLE BUILDERS
# These functions generate the CSV tables that accompany the plots.
# ==============================================================================

generate_detailed_taxon_impact_table_one_study <- function(study_data, verdicts_row) {
  study_name <- study_data$study_name
  head_prop <- if (!is.null(study_data$head_proportion)) study_data$head_proportion else 0.5
  raw_otu_table <- study_data$otu_table
  otu_tables_per_level <- study_data$otu_tables_per_level
  output_tables <- list()

  for (level in get_levels()) {
    verdict_row <- verdicts_row %>% dplyr::filter(Level == level)
    if (nrow(verdict_row) == 0) next
    safe_threshold_bp <- get_first_degradation_threshold(verdict_row)
    final_trim_bp <- if (is.na(safe_threshold_bp)) 0 else safe_threshold_bp

    # Calculate Abundances at every step
    abundance_data_list <- list()
    trim_steps <- names(otu_tables_per_level)

    for (trim_val in trim_steps) {
      current_df <- otu_tables_per_level[[trim_val]][[level]]
      if(is.null(current_df) || nrow(current_df) == 0) next
      numeric_cols <- names(current_df)[sapply(current_df, is.numeric)]
      total_rel_abun <- numeric(nrow(current_df))

      for (col in numeric_cols) {
        col_sum <- sum(current_df[[col]], na.rm = TRUE)
        if (col_sum > 0) total_rel_abun <- total_rel_abun + (current_df[[col]] / col_sum)
      }
      mean_abundance <- if (length(numeric_cols) > 0) total_rel_abun / length(numeric_cols) else 0
      abundance_data_list[[trim_val]] <- tibble::tibble(Taxon = current_df[[level]], Abundance = mean_abundance, Trim_BP = as.numeric(trim_val))
    }

    if(length(abundance_data_list) == 0) next
    full_abundance_data <- dplyr::bind_rows(abundance_data_list)
    if (nrow(full_abundance_data) == 0) next

    initial_abundances <- full_abundance_data %>% dplyr::filter(Trim_BP == 0)
    final_abundances <- full_abundance_data %>% dplyr::filter(Trim_BP == final_trim_bp)

    # Calculate Prevalence (for identifying Core Taxa)
    sample_cols <- names(raw_otu_table)[sapply(raw_otu_table, is.numeric)]
    taxonomy_map <- raw_otu_table %>%
      dplyr::select(OTU, taxonomy) %>%
      tidyr::separate(taxonomy, into = c("Domain", get_levels(), "Species"), sep = ";", fill = "right")

    for (col in get_levels()) {
      if (col %in% names(taxonomy_map)) taxonomy_map[[col]] <- sub("^[dpcofgs]__", "", taxonomy_map[[col]])
    }

    prevalence_at_level <- raw_otu_table %>%
      dplyr::left_join(taxonomy_map, by = "OTU") %>%
      dplyr::filter(!is.na(.data[[level]]) & .data[[level]] != "") %>%
      dplyr::select(all_of(level), all_of(sample_cols)) %>%
      dplyr::group_by(.data[[level]]) %>%
      dplyr::summarise(dplyr::across(all_of(sample_cols), ~sum(.x > 0) > 0), .groups = "drop") %>%
      dplyr::mutate(Prevalence_Initial = rowSums(dplyr::select(., all_of(sample_cols))) / length(sample_cols)) %>%
      dplyr::rename(Taxon = !!sym(level))

    combined_data <- initial_abundances %>%
      dplyr::select(Taxon, Abundance) %>% dplyr::rename(Abundance_Initial = Abundance) %>%
      dplyr::left_join(final_abundances %>% dplyr::select(Taxon, Abundance) %>% dplyr::rename(Abundance_Final = Abundance), by = "Taxon") %>%
      dplyr::left_join(prevalence_at_level %>% dplyr::select(Taxon, Prevalence_Initial), by = "Taxon") %>%
      dplyr::mutate(Abundance_Final = ifelse(is.na(Abundance_Final), 0, Abundance_Final),
                    Absolute_Change_Abundance = Abundance_Final - Abundance_Initial)

    # Classify Taxa for the Report
    most_abundant <- combined_data %>% dplyr::arrange(desc(Abundance_Initial)) %>% dplyr::slice_head(n = 50) %>% dplyr::mutate(Category = "Most Abundant", Rank = row_number())
    most_impacted <- combined_data %>% dplyr::arrange(desc(abs(Absolute_Change_Abundance))) %>% dplyr::slice_head(n = 50) %>% dplyr::mutate(Category = "Most Impacted", Rank = row_number())
    core_taxa <- combined_data %>% dplyr::filter(Prevalence_Initial >= 0.75) %>% dplyr::arrange(desc(Abundance_Initial)) %>% dplyr::slice_head(n = 20) %>% dplyr::mutate(Category = "Core Taxa", Rank = row_number())

    level_summary <- dplyr::bind_rows(most_abundant, most_impacted, core_taxa) %>%
      dplyr::distinct(Taxon, Category, .keep_all = TRUE) %>%
      dplyr::mutate(Study = study_name, Level = level, Head_Proportion = head_prop) %>%
      dplyr::select(Study, Level, Head_Proportion, Category, Rank, Taxon, Abundance_Initial, Abundance_Final, Absolute_Change_Abundance, Prevalence_Initial)

    output_tables[[level]] <- level_summary
  }
  dplyr::bind_rows(output_tables)
}

generate_retention_summary_table_one_study <- function(study_data, verdicts_row) {
  study_name <- study_data$study_name
  head_prop <- if (!is.null(study_data$head_proportion)) study_data$head_proportion else 0.5
  retention_summary_list <- list()

  get_counts <- function(otu_table) {
    if (is.null(otu_table) || nrow(otu_table) == 0) return(list(ASVs = 0, Reads = 0))
    total_reads <- 0
    numeric_cols <- names(otu_table)[sapply(otu_table, is.numeric)]
    for(col in numeric_cols) { total_reads <- total_reads + sum(otu_table[[col]], na.rm = TRUE) }
    return(list(ASVs = nrow(otu_table), Reads = total_reads))
  }

  for (level in get_levels()) {
    verdict_row <- verdicts_row %>% dplyr::filter(Level == level)
    if (nrow(verdict_row) == 0) next
    safe_threshold_bp <- get_first_degradation_threshold(verdict_row)
    safe_threshold_bp <- if (is.na(safe_threshold_bp)) 0 else safe_threshold_bp
    consensus_threshold_bp_raw <- study_data$target_trim_bp
    consensus_threshold_bp_lookup <- floor(consensus_threshold_bp_raw / 10) * 10

    counts_initial <- get_counts(study_data$otu_tables_per_level[['0']][[level]])
    counts_safe <- get_counts(study_data$otu_tables_per_level[[as.character(safe_threshold_bp)]][[level]])
    counts_consensus <- get_counts(study_data$otu_tables_per_level[[as.character(consensus_threshold_bp_lookup)]][[level]])

    retention_summary_list[[level]] <- tibble(
      Study = study_name,
      Level = level,
      Head_Proportion = head_prop,
      Safe_Threshold_BP = safe_threshold_bp,
      Consensus_Threshold_BP = consensus_threshold_bp_raw,
      ASVs_Initial = counts_initial$ASVs, Reads_Initial = counts_initial$Reads,
      ASVs_Retained_at_Safe = counts_safe$ASVs, Reads_Retained_at_Safe = counts_safe$Reads,
      ASVs_Retained_at_Consensus = counts_consensus$ASVs, Reads_Retained_at_Consensus = counts_consensus$Reads,
      ASV_Retention_Percent_Safe = (counts_safe$ASVs / counts_initial$ASVs) * 100,
      Read_Retention_Percent_Safe = (counts_safe$Reads / counts_initial$Reads) * 100,
      ASV_Retention_Percent_Consensus = (counts_consensus$ASVs / counts_initial$ASVs) * 100,
      Read_Retention_Percent_Consensus = (counts_consensus$Reads / counts_initial$Reads) * 100
    )
  }
  dplyr::bind_rows(retention_summary_list)
}

generate_alpha_diversity_summary_one_study <- function(study_data, verdicts_row) {
  study_name <- study_data$study_name
  head_prop <- if (!is.null(study_data$head_proportion)) study_data$head_proportion else 0.5
  alpha_div_list <- list()

  calculate_alpha <- function(otu_table) {
    if (is.null(otu_table) || nrow(otu_table) == 0) return(list(Observed = 0, Shannon = 0))
    numeric_cols <- names(otu_table)[sapply(otu_table, is.numeric)]
    if(length(numeric_cols) == 0) return(list(Observed = 0, Shannon = 0))
    vegan_otu_table <- t(otu_table[, numeric_cols, drop = FALSE])
    if (sum(vegan_otu_table) == 0) return(list(Observed = 0, Shannon = 0))
    return(list(Observed = mean(vegan::specnumber(vegan_otu_table)), Shannon = mean(vegan::diversity(vegan_otu_table, index = "shannon"))))
  }

  for (level in get_levels()) {
    verdict_row <- verdicts_row %>% dplyr::filter(Level == level)
    if (nrow(verdict_row) == 0) next
    safe_threshold_bp <- get_first_degradation_threshold(verdict_row)
    safe_threshold_bp <- if (is.na(safe_threshold_bp)) 0 else safe_threshold_bp
    consensus_threshold_bp_raw <- study_data$target_trim_bp
    consensus_threshold_bp_lookup <- floor(consensus_threshold_bp_raw / 10) * 10

    alpha_initial <- calculate_alpha(study_data$otu_tables_per_level[['0']][[level]])
    alpha_safe <- calculate_alpha(study_data$otu_tables_per_level[[as.character(safe_threshold_bp)]][[level]])
    alpha_consensus <- calculate_alpha(study_data$otu_tables_per_level[[as.character(consensus_threshold_bp_lookup)]][[level]])

    alpha_div_list[[level]] <- tibble(
      Study = study_name,
      Level = level,
      Head_Proportion = head_prop,
      Safe_Threshold_BP = safe_threshold_bp,
      Consensus_Threshold_BP = consensus_threshold_bp_raw,
      Observed_Features_Initial = alpha_initial$Observed,
      Observed_Features_Safe = alpha_safe$Observed,
      Observed_Features_Consensus = alpha_consensus$Observed,
      Shannon_Index_Initial = alpha_initial$Shannon,
      Shannon_Index_Safe = alpha_safe$Shannon,
      Shannon_Index_Consensus = alpha_consensus$Shannon
    )
  }
  dplyr::bind_rows(alpha_div_list)
}

# ------------------------------------------------------------------------------
# Function: combine_csv_streaming
# Description: Stream-combines multiple CSV files into one master CSV.
#              This uses file connections to read/write line-by-line, avoiding
#              loading the entire dataset into memory.
# ------------------------------------------------------------------------------
combine_csv_streaming <- function(pattern, output_filename, search_dir) {
  log_and_flush(paste("--- Combining files for pattern:", pattern, "---"))

  output_path <- file.path(search_dir, output_filename)
  files_to_combine <- list.files(path = search_dir, pattern = paste0("^", pattern, ".*\\.csv$"), full.names = TRUE)
  # Prevent reading the output file if it already exists from a previous run
  files_to_combine <- files_to_combine[files_to_combine != output_path]

  if (length(files_to_combine) == 0) {
    log_and_flush(paste("⚠️ No files found for pattern '", pattern, "'. Skipping summary table:", output_filename))
    return(invisible(NULL))
  }

  log_and_flush(paste("Found", length(files_to_combine), "files to combine."))

  out_conn <- file(output_path, "wt")
  tryCatch({
    # Handle First File (Include Header)
    first_file_conn <- file(files_to_combine[1], "rt")
    header <- readLines(first_file_conn, n = 1)
    writeLines(header, out_conn)
    while(length(lines <- readLines(first_file_conn, n = 10000)) > 0) {
      writeLines(lines, out_conn)
    }
    close(first_file_conn)

    # Handle Subsequent Files (Skip Header)
    if (length(files_to_combine) > 1) {
      for (f in files_to_combine[-1]) {
        in_conn <- file(f, "rt")
        readLines(in_conn, n = 1) # Skip header
        while(length(lines <- readLines(in_conn, n = 10000)) > 0) {
          writeLines(lines, out_conn)
        }
        close(in_conn)
      }
    }
  }, finally = {
    if(isOpen(out_conn)) close(out_conn)
  })

  log_and_flush(paste("✅ Successfully created combined summary:", output_filename))
}

# =========================================================================
# SECTION 4: MAIN LOGIC
# =========================================================================
main <- function() {
  log_and_flush("--- SCRIPT STARTED: Loading libraries ---")
  suppressPackageStartupMessages({
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(tidyr))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(readr))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(purrr))
    suppressPackageStartupMessages(library(tibble))
    suppressPackageStartupMessages(library(here))
    suppressPackageStartupMessages(library(changepoint))
    suppressPackageStartupMessages(library(vegan))
  })

  log_and_flush("--- Sourcing configuration and utilities ---")
  source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_config.R"))
  source(file.path(Sys.getenv("SECAT_PROJECTDIR", getwd()), "R/secat_utils.R"))

  BASE_OUTPUT_DIR <- file.path(OUTDIR)
  AGGREGATED_DATA_DIR <- file.path(BASE_OUTPUT_DIR, "aggregated_data")
  REAL_DATA_RESULTS_DIR <- file.path(BASE_OUTPUT_DIR, "real_data_results")
  SIMULATION_RESULTS_DIR <- file.path(BASE_OUTPUT_DIR, "simulation_results")

  if (!dir.exists(AGGREGATED_DATA_DIR)) dir.create(AGGREGATED_DATA_DIR, recursive = TRUE)

  # ===== STEP 1: MEMORY-EFFICIENT SIMULATION PROCESSING =====
  log_and_flush("--- STEP 1: Memory-Efficient Simulation Baseline Generation ---")
  task_dirs <- list.dirs(SIMULATION_RESULTS_DIR, recursive = FALSE, full.names = FALSE)
  task_dirs <- task_dirs[task_dirs != ""]
  log_and_flush(paste("Found", length(task_dirs), "task directories:", paste(task_dirs, collapse = ", ")))

  all_simulation_dissim_data <- list()
  if (length(task_dirs) > 0) {
    all_baseline_stats <- list()
    for (task in task_dirs) {
      task_results <- process_task_simulations_memory_safe(task, SIMULATION_RESULTS_DIR, AGGREGATED_DATA_DIR, REAL_DATA_RESULTS_DIR)
      if (!is.null(task_results) && !is.null(task_results$baseline_stats)) {
        all_baseline_stats[[task]] <- task_results$baseline_stats
        if (!is.null(task_results$full_dissim_data)) {
          all_simulation_dissim_data[[task]] <- task_results$full_dissim_data
        }
      }
      gc()
    }

    if (length(all_baseline_stats) > 0) {
      log_and_flush("Combining task-specific baselines into master file...")
      master_baseline <- bind_rows(all_baseline_stats)
      master_baseline_file <- file.path(AGGREGATED_DATA_DIR, "simulation_baseline_statistics.csv")
      write_csv(master_baseline, master_baseline_file)
      log_and_flush(paste("✅ Created master", master_baseline_file))
      combine_csv_streaming("simulation_retention_", "simulation_retention_curves.csv", AGGREGATED_DATA_DIR)
    } else {
      log_and_flush("⚠️ No valid simulation baseline data generated")
    }
  } else {
    log_and_flush("⚠️ No task directories found in simulation results")
  }

  # ===== STEP 2: PROCESS REAL DATA RESULTS =====
  log_and_flush("--- STEP 2: Processing Real Data Results ---")

  distance_cutoff       <- if (exists("DISTANCE_FROM_BASELINE_CUTOFF"))  DISTANCE_FROM_BASELINE_CUTOFF  else 0.15
  empirical_p_threshold <- if (exists("EMPIRICAL_NULL_P_THRESHOLD"))     EMPIRICAL_NULL_P_THRESHOLD     else 0.05
  penalty_method        <- if (exists("CHANGEPOINT_PENALTY_METHOD"))     toupper(CHANGEPOINT_PENALTY_METHOD) else "SIC"
  penalty_scan          <- if (exists("CHANGEPOINT_BOOTSTRAP_SCAN"))     CHANGEPOINT_BOOTSTRAP_SCAN     else seq(0.4, 2.4, by = 0.1)
  bootstrap_n           <- if (exists("CHANGEPOINT_BOOTSTRAP_N"))        CHANGEPOINT_BOOTSTRAP_N        else 50
  # New thresholds (v3.2)
  bc_ceiling            <- if (exists("BC_CEILING_THRESHOLD"))           BC_CEILING_THRESHOLD           else 0.20
  retention_floor       <- if (exists("RETENTION_FLOOR_THRESHOLD"))      RETENTION_FLOOR_THRESHOLD      else 0.70
  verdict_tolerance_bp  <- if (exists("VERDICT_CONSENSUS_TOLERANCE_BP")) VERDICT_CONSENSUS_TOLERANCE_BP else 500

  real_data_files <- list.files(path = REAL_DATA_RESULTS_DIR, pattern = "_results\\.rds$", recursive = TRUE, full.names = TRUE)
  log_and_flush(paste("Found", length(real_data_files), "real data result files."))
  if (length(real_data_files) == 0) return(invisible(NULL))

  all_verdicts <- list()
  log_and_flush("--- Starting Main Processing Loop (One Study at a Time) ---")

  # Load outlier list once at the start
  consensus_file <- file.path(OUTDIR, "intermediate/consensusregioninfo.csv")
  outlier_studies <- character(0)

  if (file.exists(consensus_file)) {
    consensus_info <- suppressMessages(readr::read_csv(consensus_file, show_col_types = FALSE))
    if (nrow(consensus_info) > 0 && !is.na(consensus_info$OutlierStudies[1]) && consensus_info$OutlierStudies[1] != "") {
      outlier_studies <- unlist(strsplit(consensus_info$OutlierStudies[1], ";"))
      outlier_studies <- trimws(outlier_studies)
      log_and_flush(paste("  -> Loaded", length(outlier_studies), "outlier studies from consensus file:", paste(outlier_studies, collapse = ", ")))
    }
  } else {
    log_and_flush("  -> No consensus file found. All studies treated as non-outliers.")
  }

  for (real_file in real_data_files) {
    current_study_data <- readRDS(real_file)
    study_name <- current_study_data$study_name
    task_id_for_sims <- if(ANALYSIS_MODE == "study") study_name else stringr::str_trim(current_study_data$primer_name)

    log_and_flush(paste("--- Processing Study:", study_name, "| Matching simulation task_id:", task_id_for_sims, "---"))

    sim_data_for_task <- tibble()
    if (task_id_for_sims %in% names(all_simulation_dissim_data)) {
      sim_data_for_task <- all_simulation_dissim_data[[task_id_for_sims]]
      log_and_flush(paste("     Found", nrow(sim_data_for_task), "simulation data points for task", task_id_for_sims))
    } else {
      log_and_flush(paste("     WARNING: No simulation data found for task", task_id_for_sims, ". Verdicts may be incomplete."))
    }

    log_and_flush("  -> Generating verdicts...")
    target_trim <- current_study_data$target_trim_bp
    real_curve <- current_study_data$dissim_data %>% arrange(Trim_BP) %>% mutate(Trim_Step = dense_rank(Trim_BP) - 1)
    sim_baseline_for_task <- sim_data_for_task %>% dplyr::filter(Trim_BP == 0)

    # --- OUTLIER LOGIC: Check once per study, before level loop ---
    is_study_outlier <- FALSE

    # Priority 1: Check if the individual study RDS has the flag (from re-runs)
    if (!is.null(current_study_data$is_outlier)) {
      is_study_outlier <- current_study_data$is_outlier
      log_and_flush(paste("     -> Outlier status from RDS:", is_study_outlier))
    }

    # Priority 2: Check the global consensus file
    if (!is_study_outlier && study_name %in% outlier_studies) {
      is_study_outlier <- TRUE
      log_and_flush(paste("     -> Identified as OUTLIER via global consensus file."))
    }

    study_verdicts <- list()
    for (level in get_levels()) {
      level_dissim <- real_curve %>% dplyr::filter(Level == level)
      cp_info    <- select_changepoint_penalty_cv(y = level_dissim$Dissimilarity, trim_bp = level_dissim$Trim_BP,
                                                  penalty_method = penalty_method, scan = penalty_scan,
                                                  nboot = bootstrap_n, min_trim_bp = 10)
      res_cutoff <- find_degradation_point(real_curve, sim_data_for_task, sim_baseline_for_task,
                                           level, task_id_for_sims, distance_cutoff, min_trim_bp = 10)
      res_null   <- find_null_model_point(real_curve, sim_data_for_task,
                                          level, task_id_for_sims, empirical_p_threshold, min_trim_bp = 10)
      res_bc     <- find_bc_ceiling_point(real_curve, level, target_trim,
                                          bc_ceiling = bc_ceiling, min_trim_bp = 10)
      res_ret    <- find_retention_floor_point(current_study_data$retention_data, level, target_trim,
                                               retention_floor = retention_floor, min_trim_bp = 10)

      # Consensus voting across all five methods
      # Consensus voting: three statistical methods only
      consensus  <- get_consensus_threshold(
        cp_bp        = cp_info$changepoint_bp,
        cut_bp       = res_cutoff$degradation_bp,
        null_bp      = res_null$nullmodel_bp,
        tolerance_bp = verdict_tolerance_bp
      )

      # Hard quality cutoffs — override verdict independently of consensus
      # If BC ceiling OR retention floor is breached, the study fails at this
      # level regardless of the consensus voting result.
      hard_fail <- (!is.na(res_bc$bc_ceiling_bp) || !is.na(res_ret$retention_floor_bp))
      if (hard_fail) {
        hard_bp  <- min(c(res_bc$bc_ceiling_bp, res_ret$retention_floor_bp), na.rm = TRUE)
        hard_max <- max(c(res_bc$bc_ceiling_bp, res_ret$retention_floor_bp), na.rm = TRUE)
        n_hard   <- sum(!is.na(c(res_bc$bc_ceiling_bp, res_ret$retention_floor_bp)))
        consensus <- list(
          status       = "HARD_FAIL",
          confirmed_bp = hard_bp,
          outer_bp     = hard_max,          # most lenient hard-gate trigger
          inner_bp     = hard_bp,           # most conservative hard-gate trigger
          n_triggered  = consensus$n_triggered + n_hard,
          n_agree      = consensus$n_agree
        )
      }

      # Pre-trim quality flag: dissimilarity already elevated at step 0
      pretrim_dissim <- real_curve %>%
        dplyr::filter(Level == level, Trim_BP == 0) %>%
        dplyr::pull(Dissimilarity)
      pretrim_flag <- length(pretrim_dissim) > 0 && !is.na(pretrim_dissim[1]) &&
                      pretrim_dissim[1] > bc_ceiling

      primer_name_for_verdict <- if(ANALYSIS_MODE == "primer") task_id_for_sims else {
        coords_file <- file.path(OUTDIR, "intermediate/study_alignment_coords.csv")
        if(file.exists(coords_file)) {
          coords_data <- read_csv(coords_file, show_col_types = FALSE)
          matched_coords <- coords_data %>% filter(study_name == !!study_name)
          if(nrow(matched_coords) > 0) matched_coords$primer_name[1] else "Unknown"
        } else "Unknown"
      }

      study_verdicts[[level]] <- tibble(
        Study                          = study_name,
        Primer                         = primer_name_for_verdict,
        Level                          = level,
        Is_Outlier                     = is_study_outlier,
        Pre_Trim_Quality_Flag          = pretrim_flag,
        Head_Proportion                = if (!is.null(current_study_data$head_proportion)) current_study_data$head_proportion else 0.5,
        # Raw method thresholds
        Threshold_Observed_Changepoint = cp_info$changepoint_bp,
        Penalty_Used_Changepoint       = cp_info$penalty_mult,
        Threshold_Observed_Cutoff      = res_cutoff$degradation_bp,
        Empirical_p_value_Cutoff       = res_cutoff$empirical_p_value,
        Threshold_Observed_NullModel   = res_null$nullmodel_bp,
        Empirical_p_value_NullModel    = res_null$empirical_p_value,
        Threshold_Observed_BCCeiling   = res_bc$bc_ceiling_bp,
        Threshold_Observed_RetFloor    = res_ret$retention_floor_bp,
        # Consensus voting result
        Consensus_Status               = consensus$status,
        Confirmed_Threshold_BP         = consensus$confirmed_bp,
        Threshold_Outer_BP             = consensus$outer_bp,   # most lenient trigger (caution band edge)
        Threshold_Inner_BP             = consensus$inner_bp,   # most conservative trigger (red band edge)
        N_Methods_Triggered            = as.integer(consensus$n_triggered),
        N_Methods_Agree                = as.integer(consensus$n_agree),
        Threshold_Required             = target_trim
      )
    }

    study_verdict_table <- dplyr::bind_rows(study_verdicts)
    all_verdicts[[study_name]] <- study_verdict_table

    log_and_flush(paste("  -> Generating and writing tables for", study_name))
    taxon_impact <- generate_detailed_taxon_impact_table_one_study(current_study_data, study_verdict_table)
    if (nrow(taxon_impact) > 0) write_csv(taxon_impact, file.path(AGGREGATED_DATA_DIR, paste0("taxon_impact_", study_name, ".csv")))

    retention_summary <- generate_retention_summary_table_one_study(current_study_data, study_verdict_table)
    if (nrow(retention_summary) > 0) write_csv(retention_summary, file.path(AGGREGATED_DATA_DIR, paste0("retention_summary_", study_name, ".csv")))

    alpha_diversity <- generate_alpha_diversity_summary_one_study(current_study_data, study_verdict_table)
    if (nrow(alpha_diversity) > 0) write_csv(alpha_diversity, file.path(AGGREGATED_DATA_DIR, paste0("alpha_diversity_", study_name, ".csv")))

    log_and_flush(paste("  -> Tables for", study_name, "written."))
    log_and_flush("  -> Cleaning memory...")

    # Explicitly remove large objects
    if (exists("current_study_data")) rm(current_study_data)
    if (exists("study_verdict_table")) rm(study_verdict_table)
    if (exists("taxon_impact")) rm(taxon_impact)
    gc()
    gc()
  }

  # ===== STEP 3: FINAL AGGREGATION =====
  log_and_flush("--- STEP 3: Final Aggregation Step ---")
  if (length(all_verdicts) > 0) {
    valid_verdicts <- Filter(function(x) !is.null(x) && nrow(x) > 0, all_verdicts)
    if (length(valid_verdicts) > 0) {
        master_verdict_table <- dplyr::bind_rows(valid_verdicts)
        write_csv(master_verdict_table, file.path(AGGREGATED_DATA_DIR, "master_verdict_table.csv"))
        log_and_flush("✅ Master verdict table saved.")
    } else {
        log_and_flush("⚠️ All studies resulted in empty verdicts. Master table not generated.")
    }
  }

  combine_csv_streaming("retention_summary_", "master_retention_summary.csv", AGGREGATED_DATA_DIR)
  combine_csv_streaming("alpha_diversity_", "master_alpha_diversity_summary.csv", AGGREGATED_DATA_DIR)
  combine_csv_streaming("taxon_impact_", "master_taxon_impact_summary.csv", AGGREGATED_DATA_DIR)

  log_and_flush("✅ All summary tables successfully generated.")
  log_and_flush("--- Aggregation Job finished ---")
}  # <-- THIS CLOSING BRACE WAS MISSING

# Run Main
tryCatch(
  expr = main(),
  error = function(e) {
    log_and_flush("\n\n--- CATASTROPHIC ERROR CAUGHT ---")
    log_and_flush(conditionMessage(e))
    log_and_flush(paste(capture.output(traceback()), collapse = "\n"))
    quit(save = "no", status = 1)
  }
)
