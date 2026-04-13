#!/usr/bin/env nextflow
// ============================================================================
// SeCAT: Sequence Consensus Amplicon Trimming
// Nextflow DSL2 Main Workflow (v4.1)
// ============================================================================

nextflow.enable.dsl = 2

include { CLEAN_DATA          } from './modules/local/clean_data'
include { STUDY_MAPPING       } from './modules/local/study_mapping'
include { COLLECT_MAPS        } from './modules/local/collect_maps'
include { PREPARE_SIMS        } from './modules/local/prepare_sims'
include { RUN_SIMULATION      } from './modules/local/run_simulation'
include { ANALYSE_REAL        } from './modules/local/analyse_real'
include { AGGREGATE           } from './modules/local/aggregate'
include { GENERATE_VERDICTS   } from './modules/local/generate_verdicts'
include { GENERATE_REPORT     } from './modules/local/generate_report'
include { MERGE_DATASETS      } from './modules/local/merge_datasets'
include { VALIDATE            } from './modules/local/validate'
include { PRIMER_MAPPING      } from './modules/local/primer_mapping'
include { GENERATE_PRIMER_DBS } from './modules/local/primer_mapping'

workflow {

    // ── Preflight checks — validate all paths before submitting any jobs ──────
    if (!params.manifest) {
        error "ERROR: --manifest is required. Set manifest in params.yaml or pass --manifest path/to/manifest.tsv"
    }
    if (!params.reference_db) {
        error "ERROR: --reference_db is required. Set reference_db in params.yaml or pass --reference_db path/to/SILVA.fasta"
    }

    def manifest_file_check = file(params.manifest)
    if (!manifest_file_check.exists()) {
        error "ERROR: Manifest file not found: ${params.manifest}"
    }

    def reference_db_check = file(params.reference_db)
    if (!reference_db_check.exists()) {
        error "ERROR: Reference database not found: ${params.reference_db}\nEnsure the SILVA aligned FASTA is accessible from the submission node."
    }

    // Validate manifest study file paths exist
    def manifest_lines = manifest_file_check.readLines()
    def manifest_headers = manifest_lines[0].split('\t') as List
    def missing_files = []
    manifest_lines[1..-1].each { line ->
        def values = line.split('\t') as List
        def row = [manifest_headers, values].transpose().collectEntries { k, v -> [(k): v] }
        ['asv_fasta_path','asv_counts_path','taxonomy_path','metadata_path'].each { col ->
            if (row.containsKey(col) && row[col] && row[col] != 'NA') {
                if (!file(row[col]).exists()) {
                    missing_files << "${row.study_name}: ${col} = ${row[col]}"
                }
            }
        }
    }
    if (missing_files) {
        error "ERROR: Missing study files referenced in manifest:\n  " + missing_files.join('\n  ')
    }

    log.info """
    ============================================================
      SeCAT v4.1 (Nextflow)
    ============================================================
      Manifest        : ${params.manifest}
      Reference DB    : ${params.reference_db}
      Analysis mode   : ${params.analysis_mode}
      Simulations/set : ${params.num_simulations}
      Trim step mode  : ${params.trim_step_mode}
      Run validation  : ${params.run_validation}
      Executor        : ${workflow.profile}
      Work dir        : ${workflow.workDir}
      Output dir      : ${params.outdir}
    ============================================================
    """.stripIndent()

    // ── Step 0: Clean data ────────────────────────────────────────────────────
    // Removes chloroplasts/mitochondria, empty samples, syncs files,
    // produces secat_manifest_clean.tsv with updated paths to cleaned files.
    CLEAN_DATA(file(params.manifest))
    clean_manifest = CLEAN_DATA.out.clean_manifest.first()

    // ── Read cleaned manifest directly in Groovy ──────────────────────────────
    // Produces a queue channel of tuples: (1-based task_id, study_name, row_map)
    // Uses the cleaned manifest output from CLEAN_DATA
    indexed_studies_ch = clean_manifest.map { mfile ->
        def mlines  = mfile.readLines()
        def headers = mlines[0].split('\t') as List
        mlines[1..-1].withIndex().collect { line, idx ->
            def values = line.split('\t') as List
            def row = [headers, values].transpose()
                          .collectEntries { k, v -> [(k): v] }
            tuple(idx + 1, row.study_name, row)
        }
    }.flatMap { it }

    // Manifest rows channel for downstream steps (trim/merge)
    manifest_rows_ch = clean_manifest.map { mfile ->
        def mlines  = mfile.readLines()
        def headers = mlines[0].split('\t') as List
        mlines[1..-1].collect { line ->
            def values = line.split('\t') as List
            def row = [headers, values].transpose()
                          .collectEntries { k, v -> [(k): v] }
            tuple(row.study_name, row)
        }
    }.flatMap { it }

    if (params.analysis_mode == 'study') {

	STUDY_MAPPING(
	    indexed_studies_ch,
	    params.reference_db,
	    clean_manifest
	)
        all_mapping_parts_ch = STUDY_MAPPING.out.mapping_part.collect()
        COLLECT_MAPS(all_mapping_parts_ch)

        collected_coords = COLLECT_MAPS.out.study_coords.first()
        collected_maps   = COLLECT_MAPS.out.collected_maps.first()

    } else if (params.analysis_mode == 'primer') {

        PRIMER_MAPPING(file(params.reference_db))

        primer_coords_ch = PRIMER_MAPPING.out.primer_coords
            .splitCsv(header: true)
            .map { row -> tuple(row.primer_name, params.reference_db) }

        GENERATE_PRIMER_DBS(primer_coords_ch)

        collected_coords = PRIMER_MAPPING.out.primer_coords
        collected_maps   = PRIMER_MAPPING.out.primer_coords

    } else {
        error "ERROR: analysis_mode must be 'study' or 'primer', got: ${params.analysis_mode}"
    }

    PREPARE_SIMS(
        collected_coords,
        params.reference_db,
        clean_manifest
    )

    ANALYSE_REAL(
        indexed_studies_ch,
        collected_coords,
        PREPARE_SIMS.out.consensus_info.first(),
        clean_manifest
    )

    sim_tasks_ch = PREPARE_SIMS.out.sim_tasks_csv.first()
        .splitCsv(header: true)
        .map { row ->
            tuple(row.task_id, row.simulation_seed, row.num_steps, row.amplicon_length)
        }

    RUN_SIMULATION(
        sim_tasks_ch,
        PREPARE_SIMS.out.sim_reference_subset.first(),
        collected_coords,
        PREPARE_SIMS.out.consensus_info
    )

    all_real_results_ch = ANALYSE_REAL.out.results_rds.collect()
    all_sim_results_ch  = RUN_SIMULATION.out.results_rds.collect()

    AGGREGATE(
        all_real_results_ch,
        all_sim_results_ch,
        collected_coords,
        PREPARE_SIMS.out.consensus_info
    )

    GENERATE_VERDICTS(AGGREGATE.out.master_verdict_table)

    GENERATE_REPORT(
        AGGREGATE.out.aggregated_dir,
        GENERATE_VERDICTS.out.verdict_data
    )
    )

    if (params.auto_trim) {
        log.warn "auto_trim=true: Trimming with KEEP-only verdicts."
        TRIM_SEQUENCES(
            GENERATE_VERDICTS.out.verdict_data,
            COLLECT_MAPS.out.study_coords,
            PREPARE_SIMS.out.consensus_info.first(),
            manifest_rows_ch.collect()
        )
        MERGE_DATASETS(
            TRIM_SEQUENCES.out.standardized_fastas.collect(),
            TRIM_SEQUENCES.out.trim_summary,
            manifest_rows_ch.collect()
        )
        if (params.run_validation) {
            VALIDATE(
                MERGE_DATASETS.out.combined_fasta,
                MERGE_DATASETS.out.combined_otu,
                MERGE_DATASETS.out.combined_tax,
                MERGE_DATASETS.out.combined_meta.ifEmpty(file('NO_META')),
                MERGE_DATASETS.out.asv_mapping.ifEmpty(file('NO_MAPPING')),
                AGGREGATE.out.aggregated_dir
            )
        }
    } else {
        log.info """
        ============================================================
        Pipeline complete. Review verdicts then run:
          Rscript R/11_select_studies.R
          nextflow run main.nf --entry STANDARDIZE -resume
        ============================================================
        """.stripIndent()
    }
}

workflow STANDARDIZE {

    selected_file  = file("${params.outdir}/aggregated_data/selected_studies_for_trim.txt",
                          checkIfExists: true)
    consensus_info = file("${params.outdir}/intermediate/consensusregioninfo.csv",
                          checkIfExists: true)
    study_coords   = file("${params.outdir}/intermediate/study_alignment_coords.csv",
                          checkIfExists: true)
    aggregated_dir = file("${params.outdir}/aggregated_data", checkIfExists: true)

    manifest_file   = file(params.manifest, checkIfExists: true)
    manifest_lines  = manifest_file.readLines()
    manifest_headers = manifest_lines[0].split('\t') as List

    manifest_rows_ch = Channel.from(
        manifest_lines[1..-1].collect { line ->
            def values = line.split('\t') as List
            def row = [manifest_headers, values].transpose()
                          .collectEntries { k, v -> [(k): v] }
            tuple(row.study_name, row)
        }
    )

    TRIM_SEQUENCES(
        selected_file,
        study_coords,
        consensus_info,
        manifest_rows_ch.collect()
    )

    MERGE_DATASETS(
        TRIM_SEQUENCES.out.standardized_fastas.collect(),
        TRIM_SEQUENCES.out.trim_summary,
        manifest_rows_ch.collect()
    )

    if (params.run_validation) {
        VALIDATE(
            MERGE_DATASETS.out.combined_fasta,
            MERGE_DATASETS.out.combined_otu,
            MERGE_DATASETS.out.combined_tax,
            MERGE_DATASETS.out.combined_meta.ifEmpty(file('NO_META')),
            MERGE_DATASETS.out.asv_mapping.ifEmpty(file('NO_MAPPING')),
            aggregated_dir
        )
    }
}

