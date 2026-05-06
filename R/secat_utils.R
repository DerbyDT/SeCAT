# ==============================================================================
# SCRIPT:   secat_utils.R
# PIPELINE: SeCAT v4.1 (Sequence Consensus Amplicon Trimming)
# PURPOSE:  Core utility library — provides all shared functions for data
#           loading, community simulation, beta-diversity analysis, VSEARCH
#           clustering, taxonomic processing, and statistical identification
#           of trimming degradation points.
#
# OVERVIEW:
#   This library is the functional backbone of the SeCAT pipeline. It is
#   sourced by every analysis script and contains functions spanning:
#   - General utilities (mode calculation, safe accessors)
#   - Null model identification (empirical p-values, consecutive significance)
#   - Community simulation (lognormal SAD, PCR bias, sequencing errors, chimeras)
#   - Data loading and validation (manifest parsing, FASTA/OTU synchronisation)
#   - Beta-diversity computation (Bray-Curtis at multiple taxonomic levels)
#   - VSEARCH integration (de novo clustering, .uc file parsing)
#   - Changepoint detection (PELT algorithm wrappers)
#   - Taxonomic retention and core taxa tracking
#
# SOURCED BY:
#   - 04_prepare_sims.R, 05_sim_worker.R, 06_analyse_real.R, 07_aggregate.R,
#     08_gen_report.R, 12_trim_sequences.R, 13_merge_datasets.R, and others
# ==============================================================================

# ==============================================================================
# SECTION 1: GENERAL PURPOSE UTILITIES
# Helper functions for basic statistical operations and safe access to
# pipeline configuration variables.
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: calculate_mode
# Purpose:  Calculates the statistical mode (most frequent value) of a vector.
#
# Parameters:
#   @param x [vector] - A vector of numbers or characters.
#
# Returns:
#   @return [any] - The most frequent value. Returns NA for empty/all-NA input.
#           Ties are broken by first occurrence (which.max behaviour).
# ------------------------------------------------------------------------------

calculate_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) { return(NA) }
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}

# ------------------------------------------------------------------------------
# Function: get_levels
# Purpose:  Safely retrieves the taxonomic levels for analysis from the global
#           environment, falling back to a sensible default (Phylum-Genus).
#
# Parameters:
#   (none) - reads TAXONOMIC_LEVELS_TO_ANALYSE from global env (secat_config.R)
#
# Returns:
#   @return [character] - Vector of taxonomic rank names to iterate over.
#
# Notes:
#   Standardises which ranks are included in diversity analyses across all
#   scripts. The fallback prevents crashes when config has not been sourced.
# ------------------------------------------------------------------------------
get_levels <- function() {
  # Check if the variable from R exists.
  if (exists("TAXONOMIC_LEVELS_TO_ANALYSE")) {
    # If it does, return its value.
    return(TAXONOMIC_LEVELS_TO_ANALYSE)
  } else {
    # If not, return a safe default to prevent crashes.
    # This also makes your code more robust.
    warning("Config variable 'TAXONOMIC_LEVELS_TO_ANALYSE' not found. Using default levels.")
    return(c("Phylum", "Class", "Order", "Family", "Genus"))
  }
}

# ==============================================================================
# SECTION 2: NULL MODEL & DEGRADATION POINT IDENTIFICATION
# Core statistical logic of SeCAT. Compares observed delta-D (change in
# beta diversity from trimming real data) against a null distribution derived
# from simulated communities, enabling data-driven detection of the trimming
# depth where community structure degrades non-randomly.
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: find_null_model_point
# Purpose:  Identifies the first trimming depth at which the real delta-D
#           becomes statistically indistinguishable from simulated delta-D
#           for a sustained window of consecutive steps.
#
# Parameters:
#   @param real_curve      [data.frame] - Dissimilarity vs trim depth for real
#          data. Requires: Level, Trim_BP, Dissimilarity.
#   @param sim_dissim_data [data.frame] - Dissimilarity vs trim depth for all
#          simulations. Requires: task_id, Level, Trim_BP, Dissimilarity,
#          simulation_id.
#   @param level           [character] - Taxonomic level (e.g., "Genus").
#   @param task_id         [character] - Study/task identifier.
#   @param empirical_p_threshold  [numeric] - P-value cutoff (default 0.05).
#   @param min_consecutive_steps  [integer] - Consecutive non-significant steps
#          needed to declare the null model point (default 3).
#   @param min_trim_bp     [integer] - Minimum trim depth to consider, avoiding
#          initial instability at short trims (default 10).
#
# Returns:
#   @return [list] - nullmodel_bp (trim depth in bp) and empirical_p_value at
#           that point. Both NA if no sustained non-significant window found.
#
# Notes:
#   Algorithm: For each step, compute delta-D (change in dissimilarity from
#   previous step) for real data, then build a null distribution of delta-D
#   from all simulations at the same interval. Empirical p-value =
#   (sum(sim_deltas >= real_delta) + 1) / (N_sims + 1), following North et al.
#   (2002). A sliding window of min_consecutive_steps where all p < threshold
#   indicates the real data has entered the noise-dominated "null" phase.
# ------------------------------------------------------------------------------

find_null_model_point <- function(real_curve, sim_dissim_data, level, task_id, empirical_p_threshold = 0.05, min_consecutive_steps = 3, min_trim_bp = 10) {
    sub_real <- real_curve %>% dplyr::filter(Level == level, Trim_BP >= 0) %>% dplyr::arrange(Trim_BP)
    if (nrow(sub_real) < 2) return(list(nullmodel_bp = NA_real_, empirical_p_value = NA_real_))

    # Calculate the change in dissimilarity at each trimming step for the real data.
    sub_real <- sub_real %>% mutate(delta_D = Dissimilarity - lag(Dissimilarity, default = Dissimilarity[1]))
    steps <- sub_real %>% dplyr::filter(Trim_BP > min_trim_bp)

    if (nrow(steps) == 0) return(list(nullmodel_bp = NA_real_, empirical_p_value = NA_real_))

    p_values <- numeric(nrow(steps))
    names(p_values) <- steps$Trim_BP

    for (i in seq_len(nrow(steps))) {
        trim_bp  <- steps$Trim_BP[i]
        real_delta_D <- steps$delta_D[i]
        # Find the immediately preceding trim step to define the interval.
        prev_trim_bp <- sub_real %>% filter(Trim_BP < trim_bp) %>% pull(Trim_BP) %>% max()

        # Filter by 'task_id' instead of 'primer_name'
        # Isolate the simulation data relevant to this specific study and trim interval.
        sim_data_for_step <- sim_dissim_data %>%
          dplyr::filter(task_id == !!task_id, Level == level, Trim_BP %in% c(prev_trim_bp, trim_bp))

        if (nrow(sim_data_for_step) == 0) { p_values[i] <- NA_real_; next }

        # For each simulation, calculate the change in dissimilarity across the interval.
        # This creates the null distribution of delta_D values.
        sim_deltas <- sim_data_for_step %>%
          group_by(simulation_id) %>%
          arrange(Trim_BP) %>%
          summarise(delta_D = Dissimilarity[2] - Dissimilarity[1], .groups = "drop") %>%
          filter(!is.na(delta_D)) %>%
          pull(delta_D)

        # Calculate the one-sided empirical p-value.
        # The +1 ensures that p is never 0, a common convention.
        # Ref: North et al. (2002) "A note on the calculation of empirical p-values".
        p_values[i] <- if(length(sim_deltas) > 0) {
          (sum(sim_deltas >= real_delta_D, na.rm = TRUE) + 1) / (length(sim_deltas) + 1)
        } else {
          NA_real_
        }
    }

    # Search for a sliding window of consecutive steps that meet the threshold.
    if (length(p_values) >= min_consecutive_steps) {
        for (i in 1:(length(p_values) - min_consecutive_steps + 1)) {
            window <- p_values[i:(i + min_consecutive_steps - 1)]
            # If all p-values in the window are non-significant, we've found the point.
            if (all(!is.na(window)) && all(window < empirical_p_threshold)) {
                return(list(nullmodel_bp = as.numeric(names(p_values)[i]), empirical_p_value = p_values[i]))
            }
        }
    }
    # If no such window is found after checking all steps, return NA.
    list(nullmodel_bp = NA_real_, empirical_p_value = NA_real_)
}

# ------------------------------------------------------------------------------
# Function: find_degradation_point
# Purpose:  Identifies the first trimming depth where delta-D exceeds a fixed
#           absolute threshold, serving as a hard backstop complementary to the
#           simulation-based null model.
#
# Parameters:
#   @param real_curve          [data.frame] - See find_null_model_point.
#   @param sim_dissim_data     [data.frame] - See find_null_model_point.
#   @param sim_dissim_baseline [data.frame] - Reserved for future comparisons.
#   @param level               [character]  - Taxonomic level to analyse.
#   @param task_id             [character]  - Study/task identifier.
#   @param distance_cutoff     [numeric]    - Absolute delta-D threshold that
#          flags degradation (default 0.15).
#   @param min_trim_bp         [integer]    - Minimum trim depth to consider.
#
# Returns:
#   @return [list] - degradation_bp and empirical_p_value (for reporting).
#           Both NA if no step exceeds the cutoff.
#
# Notes:
#   Unlike find_null_model_point (relative, simulation-based), this detects any
#   single large jump in Bray-Curtis dissimilarity -- a strong indicator of
#   systematic loss of key taxa. The empirical p-value is calculated for
#   reporting but does not drive the decision.
# ------------------------------------------------------------------------------

find_degradation_point <- function(real_curve, sim_dissim_data, sim_dissim_baseline, level, task_id, distance_cutoff = 0.15, min_trim_bp = 10) {
  sub_real <- real_curve %>% dplyr::filter(Level == level, Trim_BP >= 0) %>% dplyr::arrange(Trim_BP)
  if (nrow(sub_real) < 2) return(list(degradation_bp = NA_real_, empirical_p_value = NA_real_))

  sub_real <- sub_real %>% mutate(delta_D = Dissimilarity - lag(Dissimilarity, default = Dissimilarity[1]))
  steps <- sub_real %>% dplyr::filter(Trim_BP > min_trim_bp)

  for (i in seq_len(nrow(steps))) {
    trim_bp  <- steps$Trim_BP[i]
    distance_real <- steps$delta_D[i]
    prev_trim_bp <- sub_real %>% filter(Trim_BP < trim_bp) %>% pull(Trim_BP) %>% max()

    # Filter by 'task_id' instead of 'primer_name'
    sim_data_for_step <- sim_dissim_data %>%
      dplyr::filter(task_id == !!task_id, Level == level, Trim_BP %in% c(prev_trim_bp, trim_bp))

    if (nrow(sim_data_for_step) < 2) next

    sim_deltas <- sim_data_for_step %>%
      group_by(simulation_id) %>%
      arrange(Trim_BP) %>%
      summarise(delta_D = Dissimilarity[2] - Dissimilarity[1], .groups = "drop") %>%
      filter(!is.na(delta_D)) %>%
      pull(delta_D)

    # Calculate p-value for reporting, even though it's not used for the decision.
    empirical_p <- if(length(sim_deltas) > 0) {
      (sum(sim_deltas >= distance_real, na.rm = TRUE) + 1) / (length(sim_deltas) + 1)
    } else {
      NA_real_
    }

    # The primary condition: does the real change in dissimilarity exceed the fixed cutoff?
    if (!is.na(distance_real) && distance_real > distance_cutoff) {
      return(list(degradation_bp = trim_bp, empirical_p_value = empirical_p))
    }
  }
  # If no step exceeds the cutoff, return NA.
  list(degradation_bp = NA_real_, empirical_p_value = NA_real_)
}

# ==============================================================================
# SECTION 3: COMMUNITY SIMULATION & MANIPULATION
# Functions to generate synthetic microbial communities with realistic technical
# artifacts (PCR bias, sequencing errors, chimeras) and to process VSEARCH
# clustering outputs (.uc files) into aggregated OTU tables.
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: get_community
# Purpose:  Generates a synthetic microbial community with realistic technical
#           artifacts (PCR bias, sequencing errors, chimeras) for use as the
#           SeCAT null model.
#
# Parameters:
#   @param db_seq  [DNAStringSet] - SILVA reference database sequences.
#   @param ntaxas  [integer]      - Number of taxa to sample (default 100).
#   @param seed    [integer|NULL] - Random seed for reproducibility.
#
# Returns:
#   @return [list] with elements:
#     - sequences  [DNAStringSet] - Sequences with noise applied
#     - abundances [named numeric] - Relative abundances (sum to 1.0)
#     - table      [data.frame]   - OTU table (OTU, abundance, taxonomy)
#     - metrics    [list]         - n_taxa, shannon, evenness, n_chimeras
#
# Notes:
#   Grinder-inspired simulation pipeline (Angly et al. 2012):
#     1. Random sampling from SILVA without replacement
#     2. Lognormal SAD (mu=5, sigma=2; Preston 1948) or uniform/powerlaw
#     3. PCR bias: GC-dependent efficiency compounded over 25 cycles
#        efficiency = 1 - bias_strength * ((GC - 0.5) / 0.15)^2
#     4. Illumina errors: 0.3% substitution rate (Schirmer et al. 2015),
#        position-dependent increase toward 3' end, transition bias (70%)
#     5. Chimeras: 2% rate, random breakpoint, abundance = sqrt(p1*p2)*0.1
#        (Haas et al. 2011) -- disabled by default
#     6. Normalise to relative abundances; report Shannon H' and Pielou's J
#
#   Config variables (SIMULATION_*) control all noise model parameters.
#   Requires: Biostrings (>= 2.60.0).
# ------------------------------------------------------------------------------

get_community <- function(db_seq, ntaxas = 100, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # =================================================================
  # STEP 1: SAMPLE FROM REFERENCE DATABASE + EXTRACT TAXONOMY
  # =================================================================

  n_available <- length(db_seq)
  ntaxas <- min(ntaxas, n_available)

  if (ntaxas < 10) {
    stop(paste("Insufficient sequences in reference database:", n_available))
  }

  message(sprintf("  -> Sampling %d taxa from %d reference sequences", ntaxas, n_available))

  # Random unbiased sample
  sampled_indices <- sample(seq_along(db_seq), ntaxas, replace = FALSE)
  community_seqs <- db_seq[sampled_indices]

  # EXTRACT TAXONOMY FROM SILVA HEADERS
  # SILVA format: ">Accession.start.end Taxonomy_string"
  # Example: ">AB001234.1.1432 Bacteria;Proteobacteria;Gammaproteobacteria;..."
  silva_headers <- names(community_seqs)
  taxonomy_annotations <- sapply(silva_headers, function(header) {
    # Split on first space: "Accession Taxonomy"
    parts <- strsplit(header, " ", fixed = TRUE)[[1]]
    
    if (length(parts) >= 2) {
      # Everything after first space is taxonomy
      taxonomy <- paste(parts[-1], collapse = " ")
      
      # Clean up common SILVA formatting
      taxonomy <- gsub(";$", "", taxonomy)  # Remove trailing semicolon
      taxonomy <- trimws(taxonomy)
      
      # Ensure standard format
      if (taxonomy == "" || !grepl(";", taxonomy)) {
        return("Bacteria;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown")
      }
      
      return(taxonomy)
    } else {
      # No taxonomy in header
      return("Bacteria;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown")
    }
  })

  # Generate clean OTU names (OTU_001, OTU_002, ...)
  otu_names <- sprintf("OTU_%03d", seq_len(ntaxas))
  names(community_seqs) <- otu_names

  message(sprintf("  -> Extracted taxonomy for %d sequences", length(taxonomy_annotations)))

  # =================================================================
  # STEP 2: GENERATE INITIAL ABUNDANCES
  # =================================================================

  abundance_model <- if(exists("SIMULATION_ABUNDANCE_MODEL")) {
    SIMULATION_ABUNDANCE_MODEL
  } else {
    "lognormal"
  }

  if (abundance_model == "lognormal") {
    mu_log <- if(exists("SIMULATION_LOGNORMAL_MU")) SIMULATION_LOGNORMAL_MU else 5
    sigma_log <- if(exists("SIMULATION_LOGNORMAL_SIGMA")) SIMULATION_LOGNORMAL_SIGMA else 2
    abundances <- rlnorm(ntaxas, meanlog = mu_log, sdlog = sigma_log)

  } else if (abundance_model == "uniform") {
    abundances <- rep(1000, ntaxas)

  } else if (abundance_model == "powerlaw") {
    ranks <- seq_len(ntaxas)
    abundances <- 10000 / (ranks ^ 1.5)

  } else {
    stop("Unknown abundance model: ", abundance_model)
  }

  # Store original (pre-PCR) abundances
  original_abundances <- abundances

  # =================================================================
  # STEP 3: APPLY PCR AMPLIFICATION BIAS (Grinder Model)
  # =================================================================

  add_pcr_bias <- if(exists("SIMULATION_ADD_PCR_BIAS")) {
    SIMULATION_ADD_PCR_BIAS
  } else {
    TRUE
  }

  if (add_pcr_bias) {
    message("  -> Modeling PCR amplification bias (GC-dependent)...")

    # Calculate GC content for each sequence
    gc_content <- sapply(community_seqs, function(seq) {
      seq_char <- as.character(seq)
      seq_char <- gsub("[-.]", "", seq_char)  # Remove gaps
      bases <- strsplit(seq_char, "")[[1]]
      gc_count <- sum(bases %in% c("G", "C"))
      total_count <- sum(bases %in% c("A", "T", "G", "C"))
      if (total_count == 0) return(0.5)
      gc_count / total_count
    })

    # PCR efficiency model
    optimal_gc <- if(exists("SIMULATION_PCR_OPTIMAL_GC")) {
      SIMULATION_PCR_OPTIMAL_GC
    } else {
      0.50
    }

    bias_strength <- if(exists("SIMULATION_PCR_GC_BIAS_STRENGTH")) {
      SIMULATION_PCR_GC_BIAS_STRENGTH
    } else {
      0.65
    }

    # Quadratic penalty: sequences far from optimal GC amplify less efficiently
    # 0.15 is the scaling factor (roughly 1 SD of typical GC distributions)
    gc_deviation <- (gc_content - optimal_gc) / 0.15
    pcr_efficiency <- 1.0 - bias_strength * (gc_deviation ^ 2)
    pcr_efficiency[pcr_efficiency < 0.1] <- 0.1  # Floor at 10% efficiency

    # Simulate PCR cycles
    n_cycles <- if(exists("SIMULATION_PCR_CYCLES")) SIMULATION_PCR_CYCLES else 25

    for (cycle in 1:n_cycles) {
      abundances <- abundances * (1 + pcr_efficiency)
    }

    bias_factor <- abundances / original_abundances

    message(sprintf("     -> GC range: %.2f - %.2f", min(gc_content), max(gc_content)))
    message(sprintf("     -> PCR bias range: %.2fx - %.2fx over %d cycles",
                    min(bias_factor), max(bias_factor), n_cycles))
  }

  # =================================================================
  # STEP 4: INTRODUCE SEQUENCING ERRORS (Illumina Model)
  # =================================================================

  add_errors <- if(exists("SIMULATION_ADD_SEQUENCING_ERRORS")) {
    SIMULATION_ADD_SEQUENCING_ERRORS
  } else {
    TRUE
  }

  if (add_errors) {
    error_rate <- if(exists("SIMULATION_ERROR_RATE")) SIMULATION_ERROR_RATE else 0.003
    indel_rate <- if(exists("SIMULATION_INDEL_RATE")) SIMULATION_INDEL_RATE else 0.00003
    position_bias <- if(exists("SIMULATION_ERROR_POSITION_BIAS")) {
      SIMULATION_ERROR_POSITION_BIAS
    } else {
      TRUE
    }

    message(sprintf("  -> Introducing sequencing errors (rate: %.4f)...", error_rate))

    noisy_seqs <- lapply(seq_along(community_seqs), function(i) {
      seq <- community_seqs[[i]]
      seq_char <- strsplit(as.character(seq), "")[[1]]
      n_bases <- length(seq_char)

      # Position-dependent error rate
      if (position_bias) {
        position_multiplier <- seq(1.0, 2.0, length.out = n_bases)
      } else {
        position_multiplier <- rep(1.0, n_bases)
      }

      local_error_rates <- error_rate * position_multiplier
      has_error <- runif(n_bases) < local_error_rates
      n_errors <- sum(has_error)

      if (n_errors > 0) {
        error_positions <- which(has_error)

        for (pos in error_positions) {
          current_base <- seq_char[pos]
          if (current_base %in% c("-", ".", "N")) next

          error_type <- sample(c("sub", "ins", "del"), 1, prob = c(0.9, 0.05, 0.05))

          if (error_type == "sub") {
            if (current_base == "A") {
              seq_char[pos] <- sample(c("G", "C", "T"), 1, prob = c(0.7, 0.15, 0.15))
            } else if (current_base == "G") {
              seq_char[pos] <- sample(c("A", "C", "T"), 1, prob = c(0.7, 0.15, 0.15))
            } else if (current_base == "C") {
              seq_char[pos] <- sample(c("T", "A", "G"), 1, prob = c(0.7, 0.15, 0.15))
            } else if (current_base == "T") {
              seq_char[pos] <- sample(c("C", "A", "G"), 1, prob = c(0.7, 0.15, 0.15))
            }

          } else if (error_type == "ins") {
            if (runif(1) < 0.5) {
              insert_base <- current_base
            } else {
              insert_base <- sample(c("A", "C", "G", "T"), 1)
            }
            seq_char <- append(seq_char, insert_base, after = pos)

          } else if (error_type == "del") {
            seq_char[pos] <- ""
          }
        }
      }

      Biostrings::DNAString(paste(seq_char, collapse = ""))
    })

    community_seqs <- Biostrings::DNAStringSet(noisy_seqs)
    names(community_seqs) <- otu_names

    message(sprintf("     -> Applied errors to %d sequences", length(community_seqs)))
  }

  # =================================================================
  # STEP 5: GENERATE CHIMERAS (OPTIONAL - DISABLED BY DEFAULT)
  # =================================================================

  add_chimeras <- if(exists("SIMULATION_ADD_CHIMERAS")) SIMULATION_ADD_CHIMERAS else FALSE
  n_chimeras <- 0  # Initialize for metrics

  if (add_chimeras && length(community_seqs) >= 2) {
    chimera_rate <- if(exists("SIMULATION_CHIMERA_RATE")) SIMULATION_CHIMERA_RATE else 0.02
    n_chimeras <- max(1, ceiling(length(community_seqs) * chimera_rate))

    message(sprintf("  -> Generating %d chimeric sequences...", n_chimeras))

    chimera_seqs <- list()
    chimera_abundances <- numeric(n_chimeras)
    chimera_taxonomies <- character(n_chimeras)

    for (i in seq_len(n_chimeras)) {
      parents <- sample(seq_along(community_seqs), 2, replace = FALSE)
      parent1 <- community_seqs[[parents[1]]]
      parent2 <- community_seqs[[parents[2]]]

      # FIXED: Use nchar() for single DNAString
      len1 <- nchar(as.character(parent1))
      len2 <- nchar(as.character(parent2))
      seq_length <- min(len1, len2)
      
      breakpoint <- sample(seq(round(seq_length * 0.2), round(seq_length * 0.8)), 1)

      chimera <- Biostrings::DNAString(paste0(
        as.character(subseq(parent1, 1, breakpoint)),
        as.character(subseq(parent2, breakpoint + 1, len2))
      ))

      chimera_seqs[[i]] <- chimera
      chimera_abundances[i] <- sqrt(abundances[parents[1]] * abundances[parents[2]]) * 0.1
      
      # Chimera taxonomy
      parent1_tax <- taxonomy_annotations[parents[1]]
      tax_parts <- strsplit(parent1_tax, ";")[[1]]
      
      if (length(tax_parts) >= 2) {
        chimera_tax <- paste(
          tax_parts[1],
          tax_parts[2],
          "Chimera_Class",
          "Chimera_Order",
          "Chimera_Family",
          "Chimera_Genus",
          paste0("Chimera_", i),
          sep = ";"
        )
      } else {
        chimera_tax <- paste0("Bacteria;Chimera_Phylum;Chimera_Class;Chimera_Order;",
                             "Chimera_Family;Chimera_Genus;Chimera_", i)
      }
      
      chimera_taxonomies[i] <- chimera_tax
    }

    # Add chimeras
    chimera_set <- Biostrings::DNAStringSet(chimera_seqs)
    chimera_names <- sprintf("OTU_%03d", (length(otu_names) + 1):(length(otu_names) + n_chimeras))
    names(chimera_set) <- chimera_names

    community_seqs <- c(community_seqs, chimera_set)
    abundances <- c(abundances, chimera_abundances)
    otu_names <- c(otu_names, chimera_names)
    taxonomy_annotations <- c(taxonomy_annotations, chimera_taxonomies)
    
    message(sprintf("     -> Chimeras marked at lower taxonomic levels"))
  }

  # =================================================================
  # STEP 6: NORMALIZE AND FORMAT OUTPUT
  # =================================================================

  # Convert to relative abundances
  abundances <- abundances / sum(abundances)
  names(abundances) <- otu_names

  # Build OTU table WITH TAXONOMY
  otu_table <- data.frame(
    OTU = otu_names,
    abundance = abundances,
    taxonomy = as.character(taxonomy_annotations),
    stringsAsFactors = FALSE
  )

  # Calculate diversity metrics for QC validation
  # Shannon H' = -sum(p_i * ln(p_i)); epsilon avoids log(0)
  shannon <- -sum(abundances * log(abundances + 1e-10))
  # Pielou's J = H' / H'_max, where H'_max = ln(S)
  evenness <- shannon / log(length(abundances))

  message(sprintf("  -> Community generated: %d taxa", length(community_seqs)))
  message(sprintf("     -> Shannon diversity: %.2f | Evenness: %.2f", shannon, evenness))

  list(
    sequences = community_seqs,
    abundances = abundances,
    table = otu_table,
    metrics = list(
      n_taxa = length(community_seqs),
      shannon = shannon,
      evenness = evenness,
      n_chimeras = n_chimeras
    )
  )
}

# ------------------------------------------------------------------------------
# Function: get_otu_tab_from_uc
# Purpose:  Parses a VSEARCH .uc clustering file and aggregates the OTU table
#           by summing counts of queries that map to the same centroid.
#
# Parameters:
#   @param otu_table [data.frame] - Pre-clustering OTU table with OTU column
#          and numeric sample columns.
#   @param uc_file   [character]  - Path to the .uc clustering output.
#
# Returns:
#   @return [data.frame|NULL] - Aggregated OTU table (one row per centroid),
#           or NULL if .uc file is missing/empty/invalid.
#
# Notes:
#   After trimming, originally distinct ASVs may become identical. VSEARCH
#   clustering groups them; this function sums their counts to show how
#   taxonomic units merge. Only 'H' (hit) records are used -- these link
#   each query to its centroid. IDs are cleaned of ;size= annotations
#   before joining to ensure robust matching.
# ------------------------------------------------------------------------------

get_otu_tab_from_uc <- function(otu_table, uc_file) {
  # --- Input Validation ---
  if (!file.exists(uc_file)) {
    message(paste("    ⚠️ UC file not found:", basename(uc_file)))
    return(NULL)
  }

  if (file.info(uc_file)$size == 0) {
    message(paste("    ⚠️ UC file is empty:", basename(uc_file)))
    return(NULL)
  }

  uc <- readr::read_tsv(uc_file, col_names = FALSE, show_col_types = FALSE)

  if (!"X1" %in% names(uc)) {
    message(paste("    ⚠️ UC file has invalid format:", basename(uc_file)))
    return(NULL)
  }

  # Filter for hit records, which represent query-to-centroid mappings.
  uc_f <- uc %>% dplyr::filter(X1 == "H")

  if (nrow(uc_f) == 0) {
    message(paste("    ⚠️ No hits found in UC file:", basename(uc_file)))
    return(NULL)
  }

  # === Helper function to extract core ID ===
  # Removes size annotations (e.g., ';size=100;') and other metadata
  # commonly appended to FASTA headers, ensuring a clean join key.
  clean_id <- function(x) {
    x <- sub("\\s.*$", "", x)
    x <- sub(";.*$", "", x)
    return(x)
  }

  # Clean the UC file names (query in col X9, target/centroid in col X10)
  uc_cleaned <- uc_f %>%
    dplyr::mutate(
      query_clean = clean_id(X9),
      target_clean = clean_id(X10)
    )

  # Clean the OTU table names
  otu_cleaned <- otu_table %>%
    dplyr::mutate(OTU_clean = clean_id(OTU))

  # Create lookup: cleaned name → original name and taxonomy.
  # This allows us to restore the original centroid metadata after aggregation.
  otu_lookup <- otu_cleaned %>%
    dplyr::select(OTU_clean, OTU, taxonomy) %>%
    dplyr::distinct(OTU_clean, .keep_all = TRUE)

  # Create group mapping using cleaned names (for potential diagnostic use).
  uc_concat <- uc_cleaned %>%
    dplyr::group_by(target_clean) %>%
    dplyr::summarise(Group = stringr::str_c(query_clean, collapse = ""), .groups = 'drop')

  # Get sample column names BEFORE adding temporary columns.
  # This is robust to non-standard column names in the OTU table.
  sample_cols <- names(otu_table)[sapply(otu_table, is.numeric)]

  # --- Aggregation Logic ---
  otutable_n <- otu_cleaned %>%
    # 1. Keep only OTUs that were successfully clustered (i.e., were queries).
    dplyr::filter(OTU_clean %in% uc_cleaned$query_clean) %>%
    # 2. Join with the UC data to find the target centroid for each query OTU.
    dplyr::left_join(dplyr::select(uc_cleaned, query_clean, target_clean),
                     by = c("OTU_clean" = "query_clean")) %>%
    # 3. Group by the centroid ID.
    dplyr::group_by(target_clean) %>%
    # 4. Sum the counts for all numeric sample columns. This is the core aggregation step.
    dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), ~sum(., na.rm = TRUE)), .groups = 'drop')

  # --- Finalization ---
  # Join back the original OTU name and taxonomy for each centroid.
  otutable_n <- otutable_n %>%
    dplyr::left_join(otu_lookup, by = c("target_clean" = "OTU_clean")) %>%
    dplyr::left_join(uc_concat, by = "target_clean") %>%
    dplyr::select(-target_clean)  # Remove temporary column.

  # Diagnostic reporting.
  if (nrow(otutable_n) == 0) {
    message(paste("    ⚠️ WARNING: Zero rows returned from", basename(uc_file)))
    message(paste("    -> UC queries:", length(unique(uc_cleaned$query_clean))))
    message(paste("    -> OTU table:", nrow(otu_table)))
  } else {
    message(paste("    ✓", nrow(otutable_n), "OTUs retained from", basename(uc_file)))
  }

  return(otutable_n)
}

# ==============================================================================
# SECTION 4: DATA LOADING & VALIDATION
# Functions for reading, synchronising, and validating the primary input data
# (ASV counts tables, FASTA files, taxonomy) specified in the master manifest.
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: sync_data_for_study
# Purpose:  Loads ASV counts, FASTA sequences, and (optional) taxonomy for a
#           single study, synchronising IDs across all files.
#
# Parameters:
#   @param job_info [list/row] - Must contain asv_fasta_path, asv_counts_path;
#          optionally taxonomy_path.
#
# Returns:
#   @return [list] with elements:
#     - otu_table    [tibble]      - Counts with OTU and taxonomy columns.
#     - sequences    [BStringSet]  - Synced with otu_table row order.
#     - has_taxonomy [logical]     - TRUE if semicolon-separated taxonomy found.
#
# Notes:
#   Taxonomy is resolved in priority order: (1) external file from manifest,
#   (2) embedded column in counts table, (3) fallback to ASV IDs as taxonomy.
#   The "#OTU ID" header common in BIOM exports is cleaned automatically.
#   Fatal errors are thrown if no overlap exists between table and FASTA IDs.
# ------------------------------------------------------------------------------

sync_data_for_study <- function(job_info) {
  asv_fasta_path <- job_info$asv_fasta_path
  asv_counts_path <- job_info$asv_counts_path
  taxonomy_path <- if("taxonomy_path" %in% names(job_info) && !is.na(job_info$taxonomy_path)) {
    job_info$taxonomy_path
  } else {
    NULL
  }

  # --- File Existence Checks ---
  if (!file.exists(asv_fasta_path)) stop("FATAL: ASV FASTA file not found at: ", asv_fasta_path)
  if (!file.exists(asv_counts_path)) stop("FATAL: ASV Counts file not found at: ", asv_counts_path)

  message("  -> Loading feature table and FASTA file...")
  real_counts_raw <- readr::read_tsv(asv_counts_path, show_col_types = FALSE, skip = 1)
  real_sequences_raw <- Biostrings::readDNAStringSet(asv_fasta_path)
  message("  -> Synchronizing IDs between feature table and FASTA...")

  # Clean the feature ID column name (often starts with '#').
  feature_id_col_name <- names(real_counts_raw)[1]
  if (startsWith(feature_id_col_name, "#")) {
      new_names <- names(real_counts_raw)
      new_names[1] <- sub("^#", "", new_names[1])
      names(real_counts_raw) <- new_names
      feature_id_col_name <- names(real_counts_raw)[1]
  }

  # --- ID Synchronization ---
  ids_from_table <- real_counts_raw[[feature_id_col_name]]
  ids_from_fasta_cleaned <- names(real_sequences_raw)
  common_ids <- intersect(ids_from_table, ids_from_fasta_cleaned)
  if (length(common_ids) == 0) stop("FATAL: No common ASV IDs found between feature table and FASTA.")
  message(paste("  -> Found", length(common_ids), "common IDs."))

  # Subset both datasets to the common set of IDs.
  real_counts_synced <- real_counts_raw %>%
    dplyr::filter(.data[[feature_id_col_name]] %in% common_ids)
  real_sequences_synced <- real_sequences_raw[names(real_sequences_raw) %in% common_ids]

  message("  -> Checking for taxonomic assignments...")
  otu_table_real <- NULL
  taxonomy_data <- NULL

  # --- Taxonomy Loading (Priority 1: External File) ---
  if (!is.null(taxonomy_path) && file.exists(taxonomy_path)) {
    message("  -> Loading external taxonomy file from manifest...")
    taxonomy_data <- tryCatch({
        readr::read_tsv(taxonomy_path, show_col_types = FALSE)
    }, error = function(e) { NULL })

    if (!is.null(taxonomy_data)) {
        # Auto-detect standard ID and taxonomy column names.
        id_col <- intersect(c("Feature ID", "ASV_ID", "#OTU ID", "OTU"), names(taxonomy_data))[1]
        tax_col <- intersect(c("Taxon", "Taxonomy", "Classification", "Lineage"), names(taxonomy_data))[1]

        if (!is.na(id_col) && !is.na(tax_col)) {
            otu_table_real <- real_counts_synced %>%
                dplyr::rename(OTU = !!rlang::sym(feature_id_col_name)) %>%
                dplyr::left_join(
                  taxonomy_data %>% dplyr::select(!!rlang::sym(id_col), !!rlang::sym(tax_col)) %>% dplyr::rename(OTU = !!rlang::sym(id_col), taxonomy = !!rlang::sym(tax_col)),
                  by = "OTU"
                )
        }
    }
  }

  # --- Taxonomy Loading (Priority 2: Embedded in Counts Table) ---
  if (is.null(otu_table_real)) {
      potential_tax_cols <- names(real_counts_synced)[grepl("taxon|taxonomy|classification|lineage", names(real_counts_synced), ignore.case = TRUE)]
      if (length(potential_tax_cols) > 0) {
          message("  -> Found embedded taxonomy column in feature table.")
          taxonomy_col <- potential_tax_cols[1]
          otu_table_real <- real_counts_synced %>%
              dplyr::rename(OTU = !!rlang::sym(feature_id_col_name), taxonomy = !!rlang::sym(taxonomy_col))
      }
  }

  # --- Taxonomy Loading (Priority 3: Fallback to ASV-level) ---
  if (is.null(otu_table_real)) {
      message("  -> ⚠️  NO TAXONOMIC ASSIGNMENTS FOUND.")
      message("  -> To enable full taxonomic analysis, add a 'taxonomy_path' to your manifest.")
      message("  -> Proceeding with ASV-level analysis...")
      otu_table_real <- real_counts_synced %>%
          dplyr::rename(OTU = !!rlang::sym(feature_id_col_name)) %>%
          # Use the ASV ID itself as the "taxonomy".
          dplyr::mutate(taxonomy = OTU)
  }

  # Ensure the taxonomy column has no NAs, falling back to OTU ID if needed.
  otu_table_real <- otu_table_real %>%
    dplyr::mutate(taxonomy = ifelse(is.na(taxonomy) | taxonomy == "", OTU, taxonomy))
  # Check if the taxonomy appears to be in a standard, rank-separated format.
  taxonomic_asvs <- sum(grepl(";", otu_table_real$taxonomy), na.rm = TRUE)

  if (taxonomic_asvs > 0) {
    message("  -> ✅ Successfully loaded taxonomy. Full taxonomic analysis will be performed.")
  } else {
    message("  -> ⚠️  No semicolon-separated taxonomy found. Analysis will be at ASV-level.")
  }

  # Ensure the sequence names match the final OTU IDs in the table for perfect sync.
  names(real_sequences_synced) <- otu_table_real$OTU

  return(list(
    otu_table = otu_table_real,
    sequences = real_sequences_synced,
    has_taxonomy = taxonomic_asvs > 0
  ))
}

# ------------------------------------------------------------------------------
# Function: validate_study_data
# Purpose:  Dry-run validation of all studies in the master manifest. Checks
#           existence of required files (FASTA, counts) and optional taxonomy.
#
# Parameters:
#   (none) - reads MASTER_MANIFEST_PATH from global environment.
#
# Returns:
#   @return [data.frame] - Manifest with added columns: fasta_exists,
#           counts_exists, taxonomy_specified, taxonomy_exists,
#           ready_for_analysis. Enables programmatic filtering/reporting.
#
# Notes:
#   Run this before submitting batch jobs to catch path errors upfront,
#   avoiding mid-pipeline failures on compute clusters.
# ------------------------------------------------------------------------------

validate_study_data <- function() {
  source("R") # This seems potentially incorrect, might need clarification. Usually `source("path/to/file.R")`.
  message("--- Validating study data and taxonomy availability ---")
  message(paste("Using manifest:", MASTER_MANIFEST_PATH))

  if (!file.exists(MASTER_MANIFEST_PATH)) {
    stop("FATAL: Master manifest not found at: ", MASTER_MANIFEST_PATH)
  }

  manifest <- readr::read_tsv(MASTER_MANIFEST_PATH, show_col_types = FALSE)

  # For each study, check if essential and optional files are present on disk.
  validation_summary <- manifest %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      fasta_exists = file.exists(asv_fasta_path),
      counts_exists = file.exists(asv_counts_path),
      taxonomy_specified = "taxonomy_path" %in% names(.) && !is.na(taxonomy_path),
      taxonomy_exists = if(taxonomy_specified) file.exists(taxonomy_path) else FALSE,
      ready_for_analysis = fasta_exists && counts_exists
    ) %>%
    dplyr::ungroup()

  # --- Print a human-readable summary report to the console ---
  for (i in 1:nrow(validation_summary)) {
    study <- validation_summary[i, ]
    message(paste("Study:", study$study_name))
    message(paste("  ✓ FASTA:", if(study$fasta_exists) "Found" else "✗ Missing"))
    message(paste("  ✓ Counts:", if(study$counts_exists) "Found" else "✗ Missing"))
    message(paste("  ✓ Taxonomy:",
                  if(study$taxonomy_exists) "Found"
                  else if(study$taxonomy_specified) "✗ Missing (specified but not found)"
                  else "Not specified (will check feature table/FASTA)"))
    message("")
  }

  total_studies <- nrow(validation_summary)
  ready_studies <- sum(validation_summary$ready_for_analysis)
  taxonomy_studies <- sum(validation_summary$taxonomy_exists)

  message(paste("Summary:", ready_studies, "of", total_studies, "studies ready for analysis"))
  message(paste("Taxonomy available for:", taxonomy_studies, "studies"))

  if (ready_studies < total_studies) {
    message("⚠️  Some studies have missing required files. Please check file paths.")
  }
  return(validation_summary)
}

# ==============================================================================
# SECTION 5: CORE ANALYTICAL FUNCTIONS
# Engines for processing simulation and real-data outputs. Transforms VSEARCH
# clustering results into structured OTU tables, computes Bray-Curtis
# dissimilarity and taxonomic retention across trim steps, and applies
# changepoint detection (PELT) for threshold identification.
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: analyze_all_taxonomic_levels
# Purpose:  Iterates through all trim steps and taxonomic levels, parsing .uc
#           files into aggregated OTU tables with backfilling for missing taxa.
#
# Scientific Context:
#   To assess trimming effects at multiple resolutions (ASV through Phylum),
#   this function reconstructs the community at every step. Backfilling
#   ensures that taxa absent after trimming are retained with 0 abundance,
#   keeping matrices dimensionally consistent for valid Bray-Curtis comparison.
#
# Parameters:
#   @param otu_table [data.frame] - Original untrimmed OTU table with taxonomy.
#   @param uc_files  [list]       - List of .uc file paths keyed by trim step.
#
# Returns:
#   @return [list] - Nested: list[[trim_step]][[taxonomic_level]] = count table.
#
# Notes:
#   Algorithm: (1) parse .uc to map ASVs to centroids, (2) aggregate counts
#   at ASV and each taxonomic rank using taxonomy_map, (3) backfill missing
#   taxa with 0 counts against step "0" reference. Taxonomy prefix removal
#   (e.g., "g__") handles both QIIME2 and SILVA naming conventions.
# ------------------------------------------------------------------------------

analyze_all_taxonomic_levels <- function(otu_table, uc_files) {
  # Add "ASV" to the list of levels to be processed
  taxonomic_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV")

  has_real_taxonomy <- any(grepl(";", otu_table$taxonomy))

  if (has_real_taxonomy) {
    # Exclude "ASV" from the list of ranks to parse from the taxonomy string
    ranks_to_parse <- taxonomic_ranks[taxonomic_ranks != "ASV"]
    # Pre-calculate taxonomy breakdown for efficiency
    taxonomy_map <- otu_table %>%
      dplyr::select(OTU, taxonomy) %>%
      tidyr::separate(taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right", extra = "drop") %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(ranks_to_parse), ~stringr::str_remove(., "^[dpcofgs]__")))
  }

  otu_tables_per_level <- list()

  for (trim_val in names(uc_files)) {
    uc_entry <- uc_files[[trim_val]]

    # --- Validation of VSEARCH Output ---
    # This block of your original code for finding the current OTU table is preserved
    if (is.null(uc_entry)) {
      message(paste0("Warning: Trim level ", trim_val, " - UC entry is NULL"))
      otu_tables_per_level[[trim_val]] <- NULL
      next
    }
    has_uc_file <- ("uc_file" %in% names(uc_entry)) || ("ucfile" %in% names(uc_entry))
    if (!has_uc_file) {
      message(paste0("Warning: Trim level ", trim_val, " - UC entry missing UC file element"))
      otu_tables_per_level[[trim_val]] <- NULL
      next
    }
    uc_file_path <- if ("uc_file" %in% names(uc_entry)) {
      uc_entry$uc_file
    } else if ("ucfile" %in% names(uc_entry)) {
      uc_entry$ucfile
    } else {
      NULL
    }
    is_valid_path <- !is.null(uc_file_path) && length(uc_file_path) == 1 && !is.na(uc_file_path) && nchar(uc_file_path) > 0

    # --- Logic for Current Table Generation ---
    if (trim_val == "0") {
      message(paste0("Processing baseline (trim = 0 bp) - using UC mapping"))
      # Force baseline to use UC file like all other steps for consistency
      if (is_valid_path && file.exists(uc_file_path)) {
        result <- get_otu_tab_from_uc(otu_table, uc_file_path)
        if (!is.null(result)) {
          current_otu_table <- result
        } else {
          message("  WARNING: Baseline UC file invalid, using original table")
          current_otu_table <- otu_table
        }
      } else {
        message("  WARNING: Baseline UC file missing, using original table")
        current_otu_table <- otu_table
      }
    } else if (!is_valid_path) {
      message(paste0("Warning: Trim level ", trim_val, " bp - No valid clustering file path"))
      message("  Creating zero-abundance table (complete community loss)")
      current_otu_table <- otu_table
      sample_cols <- names(current_otu_table)[sapply(current_otu_table, is.numeric)]
      for (col in sample_cols) {
        current_otu_table[[col]] <- 0
      }
    } else if (!file.exists(uc_file_path)) {
      message(paste0("Warning: Trim level ", trim_val, " bp - UC file not found at ", uc_file_path))
      message("  Creating zero-abundance table (clustering failed or not run)")
      current_otu_table <- otu_table
      sample_cols <- names(current_otu_table)[sapply(current_otu_table, is.numeric)]
      for (col in sample_cols) {
        current_otu_table[[col]] <- 0
      }
    } else {
      result <- get_otu_tab_from_uc(otu_table, uc_file_path)
      if (is.null(result)) {
        message(paste0("Warning: Trim level ", trim_val, " bp - UC file invalid"))
        message("  Creating zero-abundance table (empty UC file)")
        current_otu_table <- otu_table
        sample_cols <- names(current_otu_table)[sapply(current_otu_table, is.numeric)]
        for (col in sample_cols) {
          current_otu_table[[col]] <- 0
        }
      } else {
        current_otu_table <- result
      }
    }
    
    # --- Aggregation by Taxonomic Level ---
    if (!is.null(current_otu_table)) {
      otu_tables_per_level[[trim_val]] <- list()

      for (level in taxonomic_ranks) {
        # Define agg_table to ensure it's available for assignment
        agg_table <- NULL

        # Special logic for the "ASV" level (finest resolution)
        if (level == "ASV") {
            agg_table <- current_otu_table %>%
                dplyr::rename(ASV = OTU) %>%
                dplyr::select(-dplyr::any_of(c("taxonomy", "Group")))

          # === EXPANDED DIAGNOSTIC BLOCK ===
          # Extensive logging for debugging data integrity at critical steps (0 and 500bp)
            if (trim_val == "0") {
                message(paste("=== BASELINE ASV TABLE ==="))
                message(paste("  Number of ASVs:", nrow(agg_table)))
                message(paste("  First 3 names:", paste(head(agg_table$ASV, 3), collapse = ", ")))

                # Check abundances
                sample_cols <- names(agg_table)[sapply(agg_table, is.numeric)]
                first_sample <- sample_cols[1]
                message(paste("  First sample column:", first_sample))
                message(paste("  Sum in", first_sample, ":", sum(agg_table[[first_sample]], na.rm=TRUE)))
                message(paste("  First 3 values:", paste(head(agg_table[[first_sample]], 3), collapse=", ")))

                total_abundance <- sum(sapply(sample_cols, function(s) sum(agg_table[[s]], na.rm=TRUE)))
                message(paste("  TOTAL ABUNDANCE across all samples:", total_abundance))
            }

            if (trim_val == "500") {
                message(paste("=== STEP 1 (500bp trim) ASV TABLE ==="))
                message(paste("  Number of ASVs:", nrow(agg_table)))
                message(paste("  First 3 names:", paste(head(agg_table$ASV, 3), collapse = ", ")))

                # Check abundances
                sample_cols <- names(agg_table)[sapply(agg_table, is.numeric)]
                first_sample <- sample_cols[1]
                message(paste("  First sample column:", first_sample))
                message(paste("  Sum in", first_sample, ":", sum(agg_table[[first_sample]], na.rm=TRUE)))
                message(paste("  First 3 values:", paste(head(agg_table[[first_sample]], 3), collapse=", ")))

                total_abundance <- sum(sapply(sample_cols, function(s) sum(agg_table[[s]], na.rm=TRUE)))
                message(paste("  TOTAL ABUNDANCE across all samples:", total_abundance))

                if ("0" %in% names(otu_tables_per_level) && "ASV" %in% names(otu_tables_per_level[["0"]])) {
                    baseline_asvs <- otu_tables_per_level[["0"]][["ASV"]][["ASV"]]
                    message(paste("  Baseline had", length(baseline_asvs), "ASVs"))
                    message(paste("  Baseline first 3:", paste(head(baseline_asvs, 3), collapse = ", ")))

                    overlap <- length(intersect(baseline_asvs, agg_table$ASV))
                    message(paste("  Overlap between baseline and step 1:", overlap, "ASVs"))

                    # Compare abundances for same sample
                    baseline_first_sample_total <- sum(otu_tables_per_level[["0"]][["ASV"]][[first_sample]], na.rm=TRUE)
                    message(paste("  Baseline", first_sample, "total:", baseline_first_sample_total))
                }
            }
            # === END DIAGNOSTIC ===

            # === BACKFILLING LOGIC (ASV Level) ===
            # If an ASV existed at Step 0 but is missing now, add it back with 0 counts.
            if (trim_val != "0" && "0" %in% names(otu_tables_per_level)) {
                if ("ASV" %in% names(otu_tables_per_level[["0"]])) {
                    original_asvs <- otu_tables_per_level[["0"]][["ASV"]][["ASV"]]
                    observed_asvs <- agg_table[["ASV"]]
                    missing_asvs <- setdiff(original_asvs, observed_asvs)

                    if (length(missing_asvs) > 0) {
                        sample_cols <- names(agg_table)[sapply(agg_table, is.numeric)]
                        missing_df <- tibble::tibble(ASV = missing_asvs)
                        for (col in sample_cols) {
                            missing_df[[col]] <- 0
                        }
                        agg_table <- dplyr::bind_rows(agg_table, missing_df)

                        # Diagnostic for backfilling
                        if (trim_val == "500") {
                            message(paste("  BACKFILLED", length(missing_asvs), "missing ASVs with zeros"))
                        }
                    }
                }
            }
            agg_table <- agg_table %>% dplyr::arrange(ASV)

        } else {
            # Original logic for all other taxonomic ranks (Phylum, Genus, etc.)
            if (has_real_taxonomy) {
              agg_table <- current_otu_table %>%
                dplyr::left_join(dplyr::select(taxonomy_map, OTU, dplyr::all_of(level)), by = "OTU") %>%
                dplyr::filter(!is.na(.data[[level]]), .data[[level]] != "") %>%
                dplyr::select(-OTU, -dplyr::any_of(c("taxonomy", "Group"))) %>%
                dplyr::group_by(.data[[level]]) %>%
                dplyr::summarise(dplyr::across(where(is.numeric), ~sum(., na.rm = TRUE)), .groups = "drop")

              # === BACKFILLING LOGIC (Taxonomic Levels) ===
              if (trim_val != "0" && "0" %in% names(otu_tables_per_level)) {
                original_taxa_at_level <- otu_tables_per_level[["0"]][[level]][[level]]
                observed_taxa <- agg_table[[level]]
                missing_taxa <- setdiff(original_taxa_at_level, observed_taxa)

                if (length(missing_taxa) > 0) {
                  sample_cols <- names(agg_table)[sapply(agg_table, is.numeric)]
                  missing_df <- tibble::tibble(!!level := missing_taxa)
                  for (col in sample_cols) {
                    missing_df[[col]] <- 0
                  }
                  agg_table <- dplyr::bind_rows(agg_table, missing_df)
                }
              }
              agg_table <- agg_table %>% dplyr::arrange(!!rlang::sym(level))

            } else {
              # Fallback for when no real taxonomy is present (treats OTU IDs as taxa)
              agg_table <- current_otu_table %>%
                dplyr::rename(!!level := OTU) %>%
                dplyr::select(-dplyr::any_of(c("taxonomy", "Group")))

              if (trim_val != "0" && "0" %in% names(otu_tables_per_level)) {
                original_taxa <- otu_tables_per_level[["0"]][[level]][[level]]
                observed_taxa <- agg_table[[level]]
                missing_taxa <- setdiff(original_taxa, observed_taxa)

                if (length(missing_taxa) > 0) {
                  sample_cols <- names(agg_table)[sapply(agg_table, is.numeric)]
                  missing_df <- tibble::tibble(!!level := missing_taxa)
                  for (col in sample_cols) {
                    missing_df[[col]] <- 0
                  }
                  agg_table <- dplyr::bind_rows(agg_table, missing_df)
                }
              }
              agg_table <- agg_table %>% dplyr::arrange(!!rlang::sym(level))
            }
        }

        # This assignment now correctly happens once per loop iteration
        otu_tables_per_level[[trim_val]][[level]] <- agg_table
      } # Closes the loop for 'level'
    } # Closes the 'if !is.null(current_otu_table)'
  } # Closes the loop for 'trim_val'

  return(otu_tables_per_level)
} # Closes the function definition

# ------------------------------------------------------------------------------
# Function: calculate_max_valid_trim
# Purpose:  Determines the deepest trim step that produced non-NULL data.
#
# Parameters:
#   @param otu_tables_per_level [list]    - Aggregated data structure.
#   @param increment            [numeric] - Step size in bp (default from config).
#
# Returns:
#   @return [numeric] - Maximum trim step number (e.g., 30 for 300bp at 10bp
#           increments). Returns 0 if no valid data exists.
# ------------------------------------------------------------------------------

calculate_max_valid_trim <- function(otu_tables_per_level, increment = NULL) {
    if (is.null(increment)) {
        increment <- if (exists("TRIM_INCREMENT")) TRIM_INCREMENT else 10
    }

    if (length(otu_tables_per_level) == 0) return(0)

    valid_trims <- names(otu_tables_per_level)[!sapply(otu_tables_per_level, is.null)]

    if (length(valid_trims) == 0) return(0)

    max_trim_bp <- max(as.numeric(valid_trims), na.rm = TRUE)
    max_trim_step <- max_trim_bp / increment

    return(max_trim_step)
}

# ------------------------------------------------------------------------------
# Function: calculate_dissimilarity
# Purpose:  Computes Bray-Curtis dissimilarity between the baseline (step 0)
#           community and the community at each successive trim step.
#
# Parameters:
#   @param otu_tables_per_level [list]    - From analyze_all_taxonomic_levels().
#   @param num_steps            [integer] - Total trim steps to analyse.
#   @param increment            [numeric] - Step size in bp (default 10).
#
# Returns:
#   @return [tibble] - Columns: Level, Trim_BP, Dissimilarity, mode.
#
# Notes:
#   This is the core metric of SeCAT. Matrices are aligned by taxon name
#   (not row order) before computing Bray-Curtis via vegan::vegdist.
#   For samples with fewer than 2 non-zero taxa, a simple Jaccard-like
#   overlap is used as fallback. Once dissimilarity reaches ~1.0, all
#   subsequent steps are set to 1.0 (total community loss).
# ------------------------------------------------------------------------------
calculate_dissimilarity <- function(otu_tables_per_level, num_steps, increment = 10) {

  dissimilarity_results <- list()

  for (level in get_levels()) {
    trim_bps <- seq(0, num_steps * increment, by = increment)
    dissim_values <- rep(NA_real_, length(trim_bps))

    original_table <- otu_tables_per_level[["0"]][[level]]

    if (is.null(original_table) || nrow(original_table) == 0) {
      message(paste("  WARNING: No OTU table for level", level, "at TrimBP=0. Skipping dissimilarity calculation."))
      next
    }

    sample_cols <- names(original_table)[sapply(original_table, is.numeric)]
    if (length(sample_cols) == 0) {
      message(paste("  WARNING: No numeric sample columns for level", level, ". Skipping dissimilarity."))
      next
    }

    original_numeric_df <- as.data.frame(original_table)[, sample_cols, drop = FALSE]
    taxa_names <- as.character(original_table[[level]])
    
    if (anyDuplicated(taxa_names)) {
      taxa_names <- make.unique(taxa_names)
    }
    rownames(original_numeric_df) <- taxa_names
    original_numeric_df[is.na(original_numeric_df)] <- 0
    
    original_filtered <- original_numeric_df[rowSums(original_numeric_df) > 0, , drop = FALSE]
    
    if (sum(original_filtered) == 0) {
      message(paste("  WARNING: All zero counts at level", level, "TrimBP=0"))
      next
    }

    dissim_values[1] <- 0
    total_loss_reached <- FALSE

    for (i in 1:num_steps) {
      if (total_loss_reached) {
        dissim_values[i + 1] <- 1.0
        next
      }

      trim_val <- i * increment
      trim_key <- as.character(trim_val)

      trimmed_table <- otu_tables_per_level[[trim_key]][[level]]

      if (is.null(trimmed_table) || nrow(trimmed_table) == 0) {
        dissim_values[i + 1] <- 1.0
        total_loss_reached <- TRUE
        next
      }

      trimmed_numeric_df <- as.data.frame(trimmed_table)[, sample_cols, drop = FALSE]
      taxa_names_trimmed <- as.character(trimmed_table[[level]])
      
      if (anyDuplicated(taxa_names_trimmed)) {
        taxa_names_trimmed <- make.unique(taxa_names_trimmed)
      }
      rownames(trimmed_numeric_df) <- taxa_names_trimmed
      trimmed_numeric_df[is.na(trimmed_numeric_df)] <- 0
      
      trimmed_filtered <- trimmed_numeric_df[rowSums(trimmed_numeric_df) > 0, , drop = FALSE]
      
      if (sum(trimmed_filtered) == 0) {
        dissim_values[i + 1] <- 1.0
        total_loss_reached <- TRUE
        next
      }

      all_taxa <- unique(c(rownames(original_filtered), rownames(trimmed_filtered)))

      original_aligned <- matrix(0, nrow = length(all_taxa), ncol = length(sample_cols), 
                                 dimnames = list(all_taxa, sample_cols))
      trimmed_aligned <- matrix(0, nrow = length(all_taxa), ncol = length(sample_cols), 
                                dimnames = list(all_taxa, sample_cols))

      taxa_in_original <- rownames(original_filtered)
      original_aligned[taxa_in_original, ] <- as.matrix(original_filtered)

      taxa_in_trimmed <- rownames(trimmed_filtered)
      trimmed_aligned[taxa_in_trimmed, ] <- as.matrix(trimmed_filtered)

      sample_dissims <- numeric(ncol(original_aligned))

      for (s in seq_len(ncol(original_aligned))) {
        orig_col <- original_aligned[, s]
        trim_col <- trimmed_aligned[, s]

        orig_nonzero <- sum(orig_col > 0)
        trim_nonzero <- sum(trim_col > 0)

        if (orig_nonzero < 2 && trim_nonzero < 2) {
          shared_taxa <- sum(orig_col > 0 & trim_col > 0)
          max_taxa <- max(orig_nonzero, trim_nonzero)
          sample_dissims[s] <- 1 - ifelse(max_taxa > 0, shared_taxa / max_taxa, 0)
        } else {
          bc_result <- tryCatch({
            vegan::vegdist(rbind(orig_col, trim_col), method = "bray")[1]
          }, error = function(e) {
            NA_real_
          })
          sample_dissims[s] <- bc_result
        }
      }

      valid_dissim <- sample_dissims[!is.na(sample_dissims)]
      
      if (length(valid_dissim) > 0) {
        dissim_values[i + 1] <- mean(valid_dissim)
      } else {
        dissim_values[i + 1] <- NA
      }
      
      if (!is.na(dissim_values[i + 1]) && dissim_values[i + 1] >= 0.9999) {
        dissim_values[i + 1] <- 1.0
        total_loss_reached <- TRUE
      }
    }

    dissimilarity_results[[level]] <- tibble::tibble(
      Level = level,
      Trim_BP = trim_bps,
      Dissimilarity = dissim_values,
      mode = "both"
    )
  }

  if (length(dissimilarity_results) == 0) return(tibble::tibble())
  
  dplyr::bind_rows(dissimilarity_results) %>%
    dplyr::arrange(Level, Trim_BP)
}

# ------------------------------------------------------------------------------
# Function: calculate_taxonomic_retention
# Purpose:  Calculates the percentage of original taxa (presence/absence)
#           retained at each trim step. Complements abundance-based Bray-Curtis.
#
# Parameters:
#   @param otu_tables_per_level [list]    - Aggregated data structure.
#   @param num_steps            [integer] - Total trim steps.
#   @param increment            [numeric] - Step size in bp (default 10).
#
# Returns:
#   @return [tibble] - Columns: Level, TrimStep, Trim_BP, Retention (%), mode.
# ------------------------------------------------------------------------------

calculate_taxonomic_retention <- function(otu_tables_per_level, num_steps, increment = 10) {
  
  all_retention <- list()
  
  for (level in get_levels()) {
    retention_values <- numeric(num_steps)
    
    # Get the initial (untrimmed) reference table
    initial_table <- otu_tables_per_level[["0"]][[level]]
    
    if (is.null(initial_table) || nrow(initial_table) == 0) {
      retention_values <- rep(0, num_steps)
    } else {
      # Identify initial taxa present (count > 0) in at least one sample
      sample_cols <- names(initial_table)[sapply(initial_table, is.numeric)]
      
      if(length(sample_cols) > 0) {
        nonzero_mask <- rowSums(initial_table[, sample_cols, drop = FALSE]) > 0
        initial_taxa <- initial_table[[level]][nonzero_mask]
        initial_taxa <- initial_taxa[!is.na(initial_taxa) & initial_taxa != ""]
        n_initial <- length(initial_taxa)
        
        for (step in 1:num_steps) {
          trim_key <- as.character(step * increment)
          trimmed_table <- otu_tables_per_level[[trim_key]][[level]]
          
          if (is.null(trimmed_table) || nrow(trimmed_table) == 0) {
            retention_values[step] <- 0
          } else {
            # Filter to taxa with nonzero abundance in any sample
            t_sample_cols <- names(trimmed_table)[sapply(trimmed_table, is.numeric)]
            if(length(t_sample_cols) > 0) {
              t_nonzero_mask <- rowSums(trimmed_table[, t_sample_cols, drop = FALSE]) > 0
              retained_taxa <- trimmed_table[[level]][t_nonzero_mask]
              retained_taxa <- retained_taxa[!is.na(retained_taxa) & retained_taxa != ""]
              
              # Count how many initial taxa are still present
              n_retained <- sum(initial_taxa %in% retained_taxa)
              retention_values[step] <- if (n_initial > 0) (n_retained / n_initial) * 100 else 0
            } else {
              retention_values[step] <- 0
            }
          }
        }
      } else {
        retention_values <- rep(0, num_steps)
      }
    }
    
    all_retention[[level]] <- tibble::tibble(
      Level = level,
      TrimStep = 1:num_steps,
      Trim_BP = (1:num_steps) * increment,
      Retention = retention_values,
      mode = "both"
    )
  }
  
  dplyr::bind_rows(all_retention)
}

# ------------------------------------------------------------------------------
# Function: calculate_changepoint_thresholds
# Purpose:  Detects tipping points in the dissimilarity curve using the PELT
#           (Pruned Exact Linear Time) changepoint algorithm.
#
# Parameters:
#   @param dissim_data [data.frame] - Output from calculate_dissimilarity().
#          Also reads CHANGEPOINT_PENALTY_METHOD and
#          CHANGEPOINT_PENALTY_MULTIPLIER from global env.
#
# Returns:
#   @return [tibble] - Columns: Level, threshold_bp.
#
# Notes:
#   Uses changepoint::cpt.meanvar with PELT to detect simultaneous shifts in
#   mean and variance. The penalty controls sensitivity: SIC (default) = BIC =
#   log(n), AIC = 2, MANUAL = user-specified. The first changepoint typically
#   marks where information loss accelerates.
# ------------------------------------------------------------------------------

calculate_changepoint_thresholds <- function(dissim_data) {
  penalty_method <- if (exists("CHANGEPOINT_PENALTY_METHOD")) toupper(CHANGEPOINT_PENALTY_METHOD) else "SIC"
  penalty_multiplier <- if (exists("CHANGEPOINT_PENALTY_MULTIPLIER")) CHANGEPOINT_PENALTY_MULTIPLIER else 1

  # Check if dissim_data is NULL or empty
  if (is.null(dissim_data) || !is.data.frame(dissim_data) || nrow(dissim_data) == 0) {
    message("  WARNING: No dissimilarity data available for changepoint analysis")
    return(tibble::tibble(
      Level = character(),
      threshold_bp = numeric()
    ))
  }

  dissim_data %>%
    dplyr::filter(!is.na(Dissimilarity)) %>%
    dplyr::group_by(Level) %>%
    dplyr::summarise(
      threshold_bp = {
        y <- Dissimilarity
        if (length(y) < 3) {
          NA_real_
        } else {
          tryCatch({
            pen_val <- switch(penalty_method,
                            "AIC" = penalty_multiplier * 2,
                            "SIC" = penalty_multiplier * log(length(y)),
                            "MANUAL" = penalty_multiplier,
                            penalty_multiplier * log(length(y)))
            # Run PELT (Pruned Exact Linear Time) algorithm
            cp_result <- changepoint::cpt.meanvar(y, method = "PELT", penalty = "Manual", pen.value = pen_val)
            cpts <- changepoint::cpts(cp_result)
            if (length(cpts) > 0 && cpts[1] <= length(y)) {
              Trim_BP[cpts[1]]
            } else {
              NA_real_
            }
          }, error = function(e) NA_real_)
        }
      },
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------
# Function: analyze_taxon_impact
# Purpose:  Tracks relative abundance of every taxon across the full trimming
#           series, enabling per-taxon trajectory plots.
#
# Parameters:
#   @param otu_tables_per_level    [list]    - Aggregated data structure.
#   @param num_steps               [integer] - Total trim steps.
#   @param increment               [numeric] - Step size in bp (default 10).
#   @param min_abundance_threshold [numeric] - Reserved for future filtering.
#
# Returns:
#   @return [list] with elements:
#     - impacted_taxa [tibble] - Long format: Taxon, Abundance, Level, Trim_BP.
#     - core_taxa     [tibble] - Placeholder for future use.
#
# Notes:
#   Abundances are normalised per-sample (column sums) then averaged across
#   samples. Missing trim steps record 0 abundance for all baseline taxa.
# ------------------------------------------------------------------------------

analyze_taxon_impact <- function(otu_tables_per_level, num_steps, increment = 10, min_abundance_threshold = 0.001) {

    all_long_format_impacts <- list()

    # Helper to correctly calculate mean relative abundance
    calculate_mean_rel_abun <- function(df, level_col_name) {
        if (is.null(df) || nrow(df) == 0) return(NULL)

        sample_cols <- names(df)[sapply(df, is.numeric)]
        if (length(sample_cols) == 0) return(NULL)

        # Normalize to relative abundance
        counts_matrix <- df[, sample_cols, drop = FALSE]
        counts_matrix[is.na(counts_matrix)] <- 0
        col_sums <- colSums(counts_matrix)
        col_sums[col_sums == 0] <- 1 # Avoid division by zero
        rel_abun_matrix <- sweep(counts_matrix, 2, col_sums, "/")

        # Create a tibble with the taxon names and their mean relative abundance
        tibble(
            Taxon = df[[level_col_name]],
            Abundance = rowMeans(rel_abun_matrix, na.rm = TRUE)
        )
    }

    for (level in get_levels()) {
        # --- Baseline (Trim Step 0) ---
        initial_table <- otu_tables_per_level[["0"]][[level]]
        if (is.null(initial_table)) next

        initial_abun <- calculate_mean_rel_abun(initial_table, level)
        if (is.null(initial_abun)) next

        # Add the baseline data for Trim_BP = 0
        all_long_format_impacts[[paste0(level, "_0")]] <- initial_abun %>%
            mutate(Level = level, Trim_BP = 0)

        # --- Loop Through All Other Trim Steps ---
        for (step in 1:num_steps) {
            trim_val <- step * increment
            trim_key <- as.character(trim_val)
            trimmed_table <- otu_tables_per_level[[trim_key]][[level]]

            # If data for this step is missing, we record zero abundance for all initial taxa
            if (is.null(trimmed_table) || nrow(trimmed_table) == 0) {
                all_long_format_impacts[[paste0(level, "_", trim_key)]] <- tibble(
                    Taxon = initial_table[[level]],
                    Abundance = 0,
                    Level = level,
                    Trim_BP = trim_val
                )
                next
            }

            trimmed_abun <- calculate_mean_rel_abun(trimmed_table, level)
            if(is.null(trimmed_abun)) next

            # Add the data for the current trim step
            all_long_format_impacts[[paste0(level, "_", trim_key)]] <- trimmed_abun %>%
                mutate(Level = level, Trim_BP = trim_val)
        }
    }

    if (length(all_long_format_impacts) == 0) {
        # Return the expected list structure, even if empty
        return(list(impacted_taxa = tibble(), core_taxa = tibble()))
    }

    # Bind all rows to create the final time-series data frame
    final_long_table <- dplyr::bind_rows(all_long_format_impacts)

    # Return the data in the list structure the report script expects
    return(list(
        impacted_taxa = final_long_table,
        core_taxa = tibble() # Keep as a placeholder for now
    ))
}

# ------------------------------------------------------------------------------
# Function: analyze_core_taxa
# Purpose:  Identifies high-prevalence ("core") taxa at baseline and tracks
#           their abundance across all trim steps for specialised reporting.
#
# Parameters:
#   @param otu_tables_per_level [list]    - Aggregated data structure.
#   @param prevalence_threshold [numeric] - Fraction of samples required for
#          "core" status (default 0.75).
#   @param top_n                [integer] - Max core taxa to report per level
#          (default 10), ranked by mean relative abundance.
#
# Returns:
#   @return [tibble|NULL] - Long format: Taxon, Abundance, Level, Trim_BP.
#           NULL if no baseline data or no taxa meet the threshold.
#
# Notes:
#   Core taxa are defined at step 0, then their trajectories are followed
#   through all trim steps. Missing taxa at a given step are recorded with
#   abundance 0 to ensure plot continuity.
# ------------------------------------------------------------------------------
analyze_core_taxa <- function(otu_tables_per_level, prevalence_threshold = 0.75, top_n = 10) {

  # Helper to get levels safely
  levels_to_process <- if(exists("get_levels")) get_levels() else c("Phylum", "Class", "Order", "Family", "Genus")

  results_list <- list()

  # 1. Identify "Core" taxa using the untrimmed (step "0") data
  baseline_data <- otu_tables_per_level[["0"]]

  if (is.null(baseline_data)) {
    warning("No baseline (step 0) data found for core taxa analysis")
    return(NULL)
  }

  for (lvl in levels_to_process) {
    if (is.null(baseline_data[[lvl]])) next

    df <- baseline_data[[lvl]]

    # Identify numeric sample columns
    sample_cols <- names(df)[sapply(df, is.numeric)]
    if (length(sample_cols) == 0) next

    # Calculate prevalence (proportion of samples where abundance > 0)
    # and mean relative abundance
    # Note: We normalize counts to relative abundance first
    counts <- df[, sample_cols, drop = FALSE]
    counts[is.na(counts)] <- 0
    col_sums <- colSums(counts)
    col_sums[col_sums == 0] <- 1
    rel_abun_mat <- sweep(counts, 2, col_sums, "/")

    prevalence <- rowMeans(counts > 0, na.rm = TRUE)
    mean_abundance <- rowMeans(rel_abun_mat, na.rm = TRUE)

    # Combine stats and filter for top N core taxa
    stats <- tibble::tibble(
      Taxon = df[[lvl]],
      Prevalence = prevalence,
      MeanAbundance = mean_abundance
    ) %>%
      dplyr::filter(!is.na(Taxon) & Taxon != "") %>%
      dplyr::filter(Prevalence >= prevalence_threshold) %>%
      dplyr::arrange(dplyr::desc(MeanAbundance)) %>%
      head(top_n)

    core_taxa_names <- stats$Taxon

    if (length(core_taxa_names) == 0) next

    # 2. Track these specific taxa across ALL trim steps
    trim_steps <- names(otu_tables_per_level)

    for (step in trim_steps) {
      step_data <- otu_tables_per_level[[step]][[lvl]]
      if (is.null(step_data)) next

      # Prepare this step's data (normalize to rel abundance)
      step_sample_cols <- names(step_data)[sapply(step_data, is.numeric)]
      if(length(step_sample_cols) == 0) next

      s_counts <- step_data[, step_sample_cols, drop = FALSE]
      s_counts[is.na(s_counts)] <- 0
      s_col_sums <- colSums(s_counts)
      s_col_sums[s_col_sums == 0] <- 1
      s_rel <- sweep(s_counts, 2, s_col_sums, "/")

      # Create lookup table for this step
      current_step_stats <- tibble::tibble(
        Taxon = step_data[[lvl]],
        Abundance = rowMeans(s_rel, na.rm = TRUE)
      ) %>%
      dplyr::filter(Taxon %in% core_taxa_names)

      # Ensure we have a row for every core taxon, even if missing (abundance 0)
      # This 'left_join' ensures continuity in plots even if a taxon drops to 0.
      full_step_stats <- tibble::tibble(Taxon = core_taxa_names) %>%
        dplyr::left_join(current_step_stats, by = "Taxon") %>%
        dplyr::mutate(
          Abundance = tidyr::replace_na(Abundance, 0),
          Level = lvl,
          Trim_BP = as.numeric(step)
        )

      results_list[[paste(lvl, step, sep = "_")]] <- full_step_stats
    }
  }

  if (length(results_list) == 0) return(NULL)

  dplyr::bind_rows(results_list)
}


# ------------------------------------------------------------------------------
# Function: cut_sequence
# Purpose:  Physically trims DNA sequences by removing a fixed number of bases
#           from the 5' (head) and/or 3' (tail) ends.
#
# Parameters:
#   @param seq  [BStringSet/DNAStringSet] - Sequences to trim.
#   @param head [integer] - Bases to remove from 5' end (default 0).
#   @param tail [integer] - Bases to remove from 3' end (default 0).
#
# Returns:
#   @return [BStringSet] - Trimmed sequences; empty strings for sequences
#           shorter than head + tail.
# ------------------------------------------------------------------------------

cut_sequence <- function(seq, head = 0, tail = 0){
  char_seq <- as.character(seq)
  new_lengths <- nchar(char_seq) - head - tail
  valid_indices <- new_lengths > 0
  x <- character(length(seq))
  x[valid_indices] <- substr(char_seq[valid_indices],
                                start = head + 1,
                                stop = nchar(char_seq[valid_indices]) - tail)
  x2 <- Biostrings::BStringSet(x)
  names(x2) <- names(seq)
  return(x2)
}

# ------------------------------------------------------------------------------
# Function: run_vsearch
# Purpose:  Wrapper for VSEARCH de novo clustering (--cluster_size).
#
# Parameters:
#   @param fasta        [character] - Input FASTA path.
#   @param out          [character] - Output prefix (produces .uc and .fasta).
#   @param vsearch_path [character] - Path to vsearch executable.
#   @param identity     [numeric]   - Clustering identity threshold (default 0.97).
#   @param strand       [character] - Strand search mode (default "both").
#
# Returns:
#   (side effect) - Writes .uc and .fasta centroid files to disk.
# ------------------------------------------------------------------------------

run_vsearch <- function(fasta, out, vsearch_path, identity = 0.97, strand = "both"){
  cmd <- sprintf('%s --cluster_size "%s" --id %s --uc "%s.uc" --strand %s --centroids "%s.fasta" --threads 1', vsearch_path, fasta, identity, out, strand, out)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
}

# ------------------------------------------------------------------------------
# Function: perform_consensus_trim
# Purpose:  Placeholder for coverage-based consensus trimming. Currently a
#           no-op that returns sequences unchanged.
#
# Parameters:
#   @param sequences       [DNAStringSet] - Input sequences.
#   @param consensus_start [integer]      - Start of consensus region.
#   @param consensus_end   [integer]      - End of consensus region.
#   @param min_coverage    [numeric]      - Minimum coverage fraction (default 0.5).
#
# Returns:
#   @return [DNAStringSet] - Sequences (currently unmodified).
# ------------------------------------------------------------------------------
perform_consensus_trim <- function(sequences, consensus_start, consensus_end, min_coverage = 0.5) {
  return(sequences)
}

# ==============================================================================
# SECTION 6: SIMULATION ENGINE
# The core trim-and-recluster loop. Progressively shortens sequences, re-maps
# them to fixed reference centroids via VSEARCH, and records which OTUs merge
# at each trim depth.
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: run_trim_analysis
# Purpose:  Core simulation engine. Progressively trims sequences in silico,
#           re-maps them to fixed reference centroids via VSEARCH, and records
#           OTU merging at each depth.
#
# Parameters:
#   @param sequences          [DNAStringSet] - Raw sequences.
#   @param vsearch_path       [character]    - Path to vsearch executable.
#   @param output_dir         [character]    - Directory for temp files/results.
#   @param num_steps          [integer]      - Number of trim steps (default 30).
#   @param mode               [character]    - "both" = proportional from each end.
#   @param increment          [numeric]      - Step size in bp (default 10).
#   @param primer_start/end   [integer|NULL] - Primer boundaries in alignment.
#   @param consensus_start/end [integer|NULL] - Consensus region boundaries.
#   @param use_alignment      [logical]      - TRUE = MSA slicing (Study Mode);
#          FALSE = physical bp cutting (Primer Mode).
#   @param aligned_sequences  [DNAStringSet] - MSA (required if use_alignment=TRUE).
#
# Returns:
#   @return [list] with elements:
#     - uc_files              [list]    - Paths to .uc files keyed by trim_bp.
#     - total_trim_to_consensus [numeric] - Total distance to consensus region.
#     - head_proportion       [numeric] - Fraction of trim applied to 5' end.
#
# Notes:
#   Study Mode: slices MSA columns to narrow the study window; robust for
#   indel-rich regions. Primer Mode: physically removes characters from
#   sequence strings. In both modes, reference centroids are established once
#   at 97% identity (step 0), then trimmed queries AND references are mapped
#   at each step. Sequences shorter than 50bp are filtered out.
#   Proportional trimming biases cuts toward whichever end is farther from
#   the consensus region.
# ------------------------------------------------------------------------------

run_trim_analysis <- function(sequences, vsearch_path, output_dir, num_steps = 30, mode = "both", increment = 10,
                              primer_start = NULL, primer_end = NULL, consensus_start = NULL, consensus_end = NULL,
                              use_alignment = FALSE, aligned_sequences = NULL) {

  temp_dir <- file.path(output_dir, "temp_vsearch")
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE)

  uc_files <- list()

  # ========================================================================
  # 1. IDENTIFY STABLE REFERENCE CENTROIDS (GROUND TRUTH)
  # ========================================================================

  reference_centroids_aligned <- NULL
  reference_centroids_physical <- NULL

  if (use_alignment && !is.null(aligned_sequences)) {
      message("  -> [Study Mode] Identifying reference centroids from aligned sequences...")

      # 1. Slice to STUDY Window (Full Coverage)
      # Isolate the region bounded by the primers (or max study window)
      study_window_aln <- Biostrings::subseq(aligned_sequences, start = primer_start, end = primer_end)

      # 2. Degap for Clustering
      # Remove gaps to perform standard identity-based clustering
      study_window_degapped <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(study_window_aln)))

      # 3. Filter empty/short sequences that are artifacts of alignment
      valid_seqs <- study_window_degapped[width(study_window_degapped) >= 50]

      if (length(valid_seqs) == 0) {
          stop("FATAL: No sequences in study region after filtering.")
      }

      message(sprintf("  -> Using %d sequences from study region [%d:%d] for centroid definition",
                      length(valid_seqs), primer_start, primer_end))

      # 4. Cluster to find Centroid IDs (97% Identity)
      seq_file_base <- file.path(temp_dir, "base_aligned_seqs.fasta")
      Biostrings::writeXStringSet(valid_seqs, seq_file_base)

      out_prefix_ref <- file.path(temp_dir, "reference")

      cmd_clust <- sprintf('%s --cluster_size "%s" --id 0.97 --centroids "%s.fasta" --uc "%s.uc" --threads 1',
                           vsearch_path, seq_file_base, out_prefix_ref, out_prefix_ref)
      system(cmd_clust, ignore.stdout = TRUE, ignore.stderr = TRUE)

      # 5. Load centroids and verify success
      centroid_file <- paste0(out_prefix_ref, ".fasta")
      if (!file.exists(centroid_file) || file.info(centroid_file)$size == 0) {
          stop("FATAL: VSEARCH clustering failed or produced no centroids.")
      }

      ref_centroids_degapped <- Biostrings::readDNAStringSet(centroid_file)

      if (length(ref_centroids_degapped) == 0) {
          stop("FATAL: VSEARCH produced empty centroid file.")
      }

      # 6. Robust name matching
      # VSEARCH may modify headers; we need to map centroids back to the original
      # aligned sequences to allow for "Alignment Mode" trimming later.
      vsearch_names <- names(ref_centroids_degapped)

      clean_vsearch_names <- sapply(vsearch_names, function(n) {
          n <- sub(";.*$", "", n)
          n <- sub("\\s.*$", "", n)
          return(n)
      })

      original_names <- names(aligned_sequences)

      clean_original_names <- sapply(original_names, function(n) {
          n <- sub(";.*$", "", n)
          n <- sub("\\s.*$", "", n)
          return(n)
      })

      matched_indices <- which(clean_original_names %in% clean_vsearch_names)

      if (length(matched_indices) == 0) {
          message("DEBUG: VSEARCH names:")
          message(paste(head(vsearch_names, 3), collapse="\n"))
          message("DEBUG: Original names:")
          message(paste(head(original_names, 3), collapse="\n"))
          stop("FATAL: Could not match any centroid names back to original sequences.")
      }

      # These aligned centroids are the parents. We will trim THEM to create the
      # reference database for each step.
      reference_centroids_aligned <- aligned_sequences[matched_indices]

      message(paste("  -> Established", length(reference_centroids_aligned), "aligned reference centroids."))

      # === CRITICAL FIX: Create proper baseline mapping (Step 0) ===
      message("  -> Creating baseline mapping (Step 0)...")

      # Define Baseline: Map all study sequences to the centroids using the full window.
      baseline_queries_aln <- Biostrings::subseq(aligned_sequences, start = primer_start, end = primer_end)
      baseline_queries <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(baseline_queries_aln)))
      baseline_queries <- baseline_queries[width(baseline_queries) >= 50]

      baseline_refs_aln <- Biostrings::subseq(reference_centroids_aligned, start = primer_start, end = primer_end)
      baseline_refs <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(baseline_refs_aln)))
      baseline_refs <- baseline_refs[width(baseline_refs) >= 50]

      # Write baseline files
      baseline_query_file <- file.path(temp_dir, "baseline_queries.fasta")
      baseline_ref_file <- file.path(temp_dir, "baseline_refs.fasta")
      Biostrings::writeXStringSet(baseline_queries, baseline_query_file)
      Biostrings::writeXStringSet(baseline_refs, baseline_ref_file)

      # Map baseline
      baseline_uc <- file.path(temp_dir, "baseline.uc")
      cmd_baseline <- sprintf('%s --usearch_global "%s" --db "%s" --id 0.97 --uc "%s" --strand both --threads 1 --maxaccepts 1 --maxrejects 32',
                             vsearch_path, baseline_query_file, baseline_ref_file, baseline_uc)
      system(cmd_baseline, ignore.stdout = TRUE, ignore.stderr = TRUE)

      uc_files[["0"]] <- list(trim_bp = 0, uc_file = baseline_uc)
      message(sprintf("  -> Baseline: %d queries, %d refs (Window: %d-%d)",
                     length(baseline_queries), length(baseline_refs), primer_start, primer_end))

  } else {
      # --- Primer Mode (Physical Trimming) ---
      message("  -> [Primer Mode] Creating reference clusters from raw sequences...")
      seq_file_ref <- file.path(temp_dir, "reference_seqs.fasta")
      Biostrings::writeXStringSet(sequences, seq_file_ref)

      out_prefix_ref <- file.path(temp_dir, "reference")
      run_vsearch(fasta = seq_file_ref, out = out_prefix_ref, vsearch_path = vsearch_path)

      reference_centroids_physical <- Biostrings::readDNAStringSet(paste0(out_prefix_ref, ".fasta"))
      uc_files[["0"]] <- list(trim_bp = 0, uc_file = paste0(out_prefix_ref, ".uc"))
  }


  # ========================================================================
  # 2. CALCULATE PROPORTIONS (Directional Trimming Logic)
  # ========================================================================
  # Determine if we should trim equally from both ends, or bias towards one end
  # based on where the consensus region sits relative to the primers.
  head_proportion <- 0.5
  total_dist_to_cons <- 0

  if (!is.null(primer_start) && !is.null(consensus_start) && !is.null(primer_end) && !is.null(consensus_end)) {
      d_h <- max(0, consensus_start - primer_start)
      d_t <- max(0, primer_end - consensus_end)

      tot <- d_h + d_t
      if (tot > 0) {
          head_proportion <- d_h / tot
          message(sprintf("  -> Proportional Logic: Head=%.2f (Dist %d), Tail=%.2f (Dist %d)",
                          head_proportion, d_h, 1.0-head_proportion, d_t))
      }
      total_dist_to_cons <- tot
  }

  # ========================================================================
  # 3. TRIM LOOP
  # ========================================================================
  mode_str <- if(use_alignment) "Alignment Slicing" else "Physical BP"
  message(sprintf("  -> Starting Trim Loop (%d steps, Mode: %s)...", num_steps, mode_str))

  MIN_SEQ_LENGTH <- 50  # Minimum length for reliable mapping

  for (i in 1:num_steps) {
      trim_bp <- i * increment

      # Calculate cut lengths for this step
      curr_head_cut <- round(trim_bp * head_proportion)
      curr_tail_cut <- trim_bp - curr_head_cut

      trimmed_queries <- NULL
      trimmed_refs <- NULL

      if (use_alignment && !is.null(reference_centroids_aligned)) {
          # --- ALIGNMENT MODE ---
          # Shrink the window by moving start/end coordinates inward
          w_start <- primer_start + curr_head_cut
          w_end   <- primer_end   - curr_tail_cut

          if (w_start >= w_end) {
              message("  -> Window closed. Stopping."); break
          }

          # 1. Slice & Degap QUERIES (The population being simulated)
          q_slice <- Biostrings::subseq(aligned_sequences, start = w_start, end = w_end)
          trimmed_queries <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(q_slice)))

          original_query_count <- length(trimmed_queries)
          trimmed_queries <- trimmed_queries[width(trimmed_queries) >= MIN_SEQ_LENGTH]

          if (length(trimmed_queries) < original_query_count * 0.1) {
              message(sprintf("  -> WARNING: Step %d - Only %d/%d sequences have >=%dbp. Stopping.",
                             i, length(trimmed_queries), original_query_count, MIN_SEQ_LENGTH))
              break
          }

          # 2. Slice & Degap REFERENCES (The target centroids)
          # We must also trim the references so we are mapping "short query" to "short ref"
          r_slice <- Biostrings::subseq(reference_centroids_aligned, start = w_start, end = w_end)
          trimmed_refs <- Biostrings::DNAStringSet(gsub("[-.]", "", as.character(r_slice)))
          trimmed_refs <- trimmed_refs[width(trimmed_refs) >= MIN_SEQ_LENGTH]

          if (length(trimmed_refs) == 0) {
              message(sprintf("  -> Step %d - No reference centroids in window. Stopping.", i))
              break
          }

          # Diagnostic message
          message(sprintf("  -> Step %d: %d queries, %d refs (Window: %d-%d, Trim: %dbp)",
                         i, length(trimmed_queries), length(trimmed_refs), w_start, w_end, trim_bp))

      } else {
          # --- PHYSICAL MODE ---
          # Just cut characters from the string
          trimmed_queries <- cut_sequence(sequences, head = curr_head_cut, tail = curr_tail_cut)
          trimmed_refs <- cut_sequence(reference_centroids_physical, head = curr_head_cut, tail = curr_tail_cut)
      }

      if (length(trimmed_queries) == 0 || length(trimmed_refs) == 0) {
          message(paste("  -> Sequences lost at step", i)); break
      }

      # Write Temp Files
      current_fasta <- file.path(temp_dir, paste0("trim_", trim_bp, ".fasta"))
      Biostrings::writeXStringSet(trimmed_queries, current_fasta)

      current_ref_path <- file.path(temp_dir, paste0("ref_trim_", trim_bp, ".fasta"))
      Biostrings::writeXStringSet(trimmed_refs, current_ref_path)

      # Map Queries -> References
      # If a query maps to its original centroid -> Success (information retained)
      # If a query maps to a different centroid -> Ambiguity/Loss of resolution
      out_uc <- file.path(temp_dir, paste0("trim_", trim_bp))
      uc_file <- paste0(out_uc, ".uc")

      cmd_map <- sprintf('%s --usearch_global "%s" --db "%s" --id 0.97 --uc "%s" --strand both --threads 1 --maxaccepts 1 --maxrejects 32',
                         vsearch_path, current_fasta, current_ref_path, uc_file)

      system(cmd_map, ignore.stdout = TRUE, ignore.stderr = TRUE)

      uc_files[[as.character(trim_bp)]] <- list(trim_bp = trim_bp, uc_file = uc_file)
  }

  return(list(
      uc_files = uc_files,
      total_trim_to_consensus = total_dist_to_cons,
      head_proportion = head_proportion
  ))
}

# =====================================================================
