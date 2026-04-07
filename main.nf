#!/usr/bin/env nextflow
// ============================================================================
// SeCAT: Sequence Consensus Amplicon Trimming
// Nextflow DSL2 Main Workflow (v4.0)
//
// USAGE:
//   # Run full pipeline (stages 1-6):
//   nextflow run main.nf -profile sge -params-file params.yaml
//
//   # Resume after failure:
//   nextflow run main.nf -profile sge -params-file params.yaml -resume
//
//   # Run stage 7 (trim + merge) after interactive study selection:
//   nextflow run main.nf -profile sge -params-file params.yaml --entry STANDARDIZE -resume
//
// PIPELINE FLOW:
//   STUDY_MAPPING (parallel x N studies)
//       └─> COLLECT_MAPS
//                ├─> PREPARE_SIMS
//                │       ├─> RUN_SIMULATION (parallel x 1100 tasks, max 50 concurrent)
//                │       └─> ANALYSE_REAL   (parallel x N studies)
//                │                 └─> AGGREGATE (waits for all sim + real results)
//                │                           ├─> GENERATE_VERDICTS
//                │                           ├─> GENERATE_REPORT
//                │                           └─> GENERATE_INDEX
//                │
//                └─> [Manual: Rscript R/11_select_studies.R]
//                        └─> nextflow run main.nf --entry STANDARDIZE -resume
//                                  ├─> TRIM_SEQUENCES
//                                  └─> MERGE_DATASETS
// ============================================================================

nextflow.enable.dsl = 2

// ----------------------------------------------------------------------------
// Import process modules
// ----------------------------------------------------------------------------
include { STUDY_MAPPING     } from './modules/local/study_mapping'
include { COLLECT_MAPS      } from './modules/local/collect_maps'
include { PREPARE_SIMS      } from './modules/local/prepare_sims'
include { RUN_SIMULATION    } from './modules/local/run_simulation'
include { ANALYSE_REAL      } from './modules/local/analyse_real'
include { AGGREGATE         } from './modules/local/aggregate'
include { GENERATE_VERDICTS } from './modules/local/generate_verdicts'
include { GENERATE_REPORT   } from './modules/local/generate_report'
include { GENERATE_INDEX    } from './modules/local/generate_index'
include { TRIM_SEQUENCES    } from './modules/local/trim_sequences'
include { MERGE_DATASETS    } from './modules/local/merge_datasets'

// ============================================================================
// MAIN WORKFLOW — Stages 1 through 6 (+ optional Stage 7 via auto_trim)
// ============================================================================

workflow {

    // -------------------------------------------------------------------------
    // Input validation
    // -------------------------------------------------------------------------
    if (!params.manifest) {
        error "ERROR: --manifest is required. Set it in params.yaml or pass --manifest secat_manifest.tsv"
    }
    if (!params.reference_db) {
        error "ERROR: --reference_db is required. Set it in params.yaml or pass --reference_db /path/to/SILVA.fasta"
    }

    log.info """
    ============================================================
      SeCAT v4.0 (Nextflow)
    ============================================================
      Manifest        : ${params.manifest}
      Reference DB    : ${params.reference_db}
      Simulations/set : ${params.num_simulations}
      Trim step mode  : ${params.trim_step_mode}
      Executor        : ${workflow.profile}
      Work dir        : ${workflow.workDir}
      Output dir      : ${params.outdir}
    ============================================================
    """.stripIndent()

    // -------------------------------------------------------------------------
    // Build manifest channel
    // Produces one element per study row: tuple(study_name, row_map)
    // -------------------------------------------------------------------------
    manifest_ch = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.study_name, row) }

    // Enumerate studies (1-based index) for R scripts that use SGE_TASK_ID.
    // .collect() gathers all rows first, then flatMap re-emits them with an
    // index — guaranteeing stable 1-based numbering regardless of emit order.
    // Produces: tuple(task_id, study_name, row_map)
    indexed_studies_ch = manifest_ch
        .collect()
        .flatMap { rows ->
            rows.withIndex().collect { item, idx ->
                tuple(idx + 1, item[0], item[1])
            }
        }

    // -------------------------------------------------------------------------
    // Stage 2: Study mapping — parallel, one task per study
    // Aligns each study's ASVs to SILVA via DECIPHER, produces:
    //   - mapping_part_N.csv           (collected by COLLECT_MAPS)
    //   - {study}_coords.csv           (per-ASV coordinates)
    //   - {study}_aligned.fasta        (aligned sequences, needed by TRIM_SEQUENCES)
    // -------------------------------------------------------------------------
    STUDY_MAPPING(
        indexed_studies_ch,
        file(params.reference_db)
    )

    // -------------------------------------------------------------------------
    // Collect all mapping parts before proceeding.
    // .collect() waits for every STUDY_MAPPING task to finish — equivalent to
    // the old SGE -hold_jid on the study mapping array job.
    // COLLECT_MAPS aggregates the per-study CSV parts into study_alignment_coords.csv
    // -------------------------------------------------------------------------
    COLLECT_MAPS(
        STUDY_MAPPING.out.mapping_part.collect()
    )

    // -------------------------------------------------------------------------
    // Collect aligned FASTAs from STUDY_MAPPING output channel.
    // Storing them here as a channel value means Nextflow tracks them as
    // explicit named inputs to TRIM_SEQUENCES rather than relying on
    // publishDir filesystem paths (which can fail with symlinks on some clusters).
    // -------------------------------------------------------------------------
    aligned_fastas_ch = STUDY_MAPPING.out.aligned_fasta.collect()

    // -------------------------------------------------------------------------
    // Stage 3: Prepare simulation tasks
    // Reads study coords, samples SILVA reference subset, writes:
    //   - simulation_tasks.csv              (N studies x num_simulations rows)
    //   - simulation_reference_subset.fasta (10k SILVA sequences for simulations)
    //   - consensusregioninfo.csv           (global consensus region boundaries)
    // -------------------------------------------------------------------------
    PREPARE_SIMS(
        COLLECT_MAPS.out.study_coords,
        file(params.reference_db)
    )

    // -------------------------------------------------------------------------
    // Stage 4: Real data analysis — parallel, one task per study
    // Runs concurrently with Stage 5 (both depend only on PREPARE_SIMS output).
    // Produces {study_name}_results.rds per study.
    // -------------------------------------------------------------------------
    ANALYSE_REAL(
        indexed_studies_ch,
        COLLECT_MAPS.out.study_coords,
        PREPARE_SIMS.out.consensus_info
    )

    // -------------------------------------------------------------------------
    // Stage 5: Simulation workers — fan-out from task channel
    //
    // splitCsv() converts simulation_tasks.csv (1100 rows) into a channel
    // emitting one tuple per row. Each tuple becomes one RUN_SIMULATION task,
    // running in parallel automatically — no array job submission needed.
    //
    // maxForks = 50 in nextflow.config throttles concurrency, equivalent to
    // the old SGE -tc 50 flag on the simulation worker array.
    // -------------------------------------------------------------------------
    sim_tasks_ch = PREPARE_SIMS.out.sim_tasks_csv
        .splitCsv(header: true)
        .map { row ->
            tuple(row.task_id, row.simulation_seed, row.num_steps, row.amplicon_length)
        }

    RUN_SIMULATION(
        sim_tasks_ch,
        PREPARE_SIMS.out.sim_reference_subset,
        COLLECT_MAPS.out.study_coords,
        PREPARE_SIMS.out.consensus_info
    )

    // -------------------------------------------------------------------------
    // Stage 6: Aggregation
    //
    // .collect() on both output channels waits for ALL real data analyses and
    // ALL simulation tasks to complete before AGGREGATE runs — equivalent to
    // the old: qsub -hold_jid "$REAL_DATA_JOB_ID,$SIM_JOB_ID" 07_aggregate.sge
    // -------------------------------------------------------------------------
    AGGREGATE(
        ANALYSE_REAL.out.results_rds.collect(),
        RUN_SIMULATION.out.results_rds.collect(),
        COLLECT_MAPS.out.study_coords,
        PREPARE_SIMS.out.consensus_info
    )

    // -------------------------------------------------------------------------
    // Stage 6b: Verdicts, reports, and HTML index
    // All three run in parallel once AGGREGATE finishes — independent of each other.
    // -------------------------------------------------------------------------
    GENERATE_VERDICTS(
        AGGREGATE.out.master_verdict_table
    )

    GENERATE_REPORT(
        AGGREGATE.out.aggregated_dir,
        GENERATE_VERDICTS.out.verdict_data
    )

    GENERATE_INDEX(
        AGGREGATE.out.aggregated_dir,
        GENERATE_VERDICTS.out.verdict_data
    )

    // -------------------------------------------------------------------------
    // Stage 7 — two paths:
    //
    // DEFAULT (auto_trim = false):
    //   Pipeline stops here. User reviews verdicts, runs 11_select_studies.R
    //   interactively, then re-launches with --entry STANDARDIZE.
    //
    // auto_trim = true:
    //   Skips human review and immediately trims all KEEP-verdict studies.
    //   Only appropriate if you have already validated the verdict logic and
    //   do not need to manually exclude any borderline studies.
    // -------------------------------------------------------------------------
    if (params.auto_trim) {

        log.warn "auto_trim = true: proceeding without human review of verdicts. " +
                 "Check ${params.outdir}/aggregated_data/verdict_data_all_levels.csv first."

        TRIM_SEQUENCES(
            GENERATE_VERDICTS.out.verdict_data,
            COLLECT_MAPS.out.study_coords,
            PREPARE_SIMS.out.consensus_info,
            manifest_ch.collect(),
            aligned_fastas_ch
        )

        MERGE_DATASETS(
            TRIM_SEQUENCES.out.standardized_fastas.collect(),
            TRIM_SEQUENCES.out.trim_summary,
            manifest_ch.collect()
        )

    } else {

        // Emit next-steps guidance to the Nextflow log when the verdict file lands
        GENERATE_VERDICTS.out.verdict_data.view { f ->
            """
            ============================================================
            SeCAT pipeline complete (Stages 1-6).

            Verdicts written to:
              ${params.outdir}/aggregated_data/verdict_data_all_levels.csv

            NEXT STEPS:
              1. Review reports:
                   ${params.outdir}/reports/

              2. Run interactive study selector (on login node):
                   Rscript R/11_select_studies.R

              3. Submit trim + merge:
                   nextflow run main.nf \\
                     -profile ${workflow.profile} \\
                     -params-file params.yaml \\
                     --entry STANDARDIZE \\
                     -resume
            ============================================================
            """.stripIndent()
        }
    }
}

// ============================================================================
// STANDARDIZE entry point — Stage 7
//
// Run this after Rscript R/11_select_studies.R has written:
//   output/aggregated_data/selected_studies_for_trim.txt
//
// USAGE:
//   nextflow run main.nf \
//     --entry STANDARDIZE \
//     -profile sge \
//     -params-file params.yaml \
//     -resume
//
// -resume is important: it means Nextflow reads aligned FASTAs from the
// cached STUDY_MAPPING outputs rather than re-running alignment.
// ============================================================================

workflow STANDARDIZE {

    // Verify required files from the main pipeline run exist before proceeding
    selected_file  = file("${params.outdir}/aggregated_data/selected_studies_for_trim.txt",
                          checkIfExists: true)
    consensus_info = file("${params.outdir}/intermediate/consensusregioninfo.csv",
                          checkIfExists: true)
    study_coords   = file("${params.outdir}/intermediate/study_alignment_coords.csv",
                          checkIfExists: true)

    // Collect aligned FASTAs from the published output directory.
    // Written by STUDY_MAPPING during the main workflow run.
    // checkIfExists: true gives a clear error if study mapping never completed.
    aligned_fastas_ch = Channel
        .fromPath("${params.outdir}/intermediate/aligned_fastas/*_aligned.fasta",
                  checkIfExists: true)
        .collect()

    manifest_ch = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.study_name, row) }

    TRIM_SEQUENCES(
        selected_file,
        study_coords,
        consensus_info,
        manifest_ch.collect(),
        aligned_fastas_ch
    )

    MERGE_DATASETS(
        TRIM_SEQUENCES.out.standardized_fastas.collect(),
        TRIM_SEQUENCES.out.trim_summary,
        manifest_ch.collect()
    )
}
