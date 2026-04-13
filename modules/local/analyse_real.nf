process ANALYSE_REAL {
    tag "${study_name}"
    label 'mem_64g'
    publishDir "${params.outdir}/real_data_results/${study_name}", mode: 'copy', pattern: "*_results.rds"

    input:
    tuple val(task_id), val(study_name), val(manifest_row)
    path study_coords
    path consensus_info
    path clean_manifest

    output:
    path "${study_name}_results.rds", emit: results_rds

    script:
    """
    export SGE_TASK_ID=${task_id}
    export SECAT_MANIFEST="\$(realpath ${clean_manifest})"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_VSEARCH_PATH="vsearch"  # binary is on PATH inside container
    export SECAT_TRIM_STEP_MODE="${params.trim_step_mode}"
    export SECAT_TRIM_INCREMENT="${params.trim_increment}"
    export SECAT_DEFAULT_MAX_TRIM_STEPS="${params.default_max_trim_steps}"
    export SECAT_CONSENSUS_BUFFER_STEPS="${params.consensus_buffer_steps}"
    export SECAT_MAX_ABSOLUTE_TRIM_STEPS="${params.max_absolute_trim_steps}"
    export SECAT_MIN_RELATIVE_ABUNDANCE="${params.min_relative_abundance}"
    export SECAT_MIN_TAXA_FOR_BRAY="${params.min_taxa_for_bray}"
    export SECAT_MIN_MEDIAN_READ_DEPTH="${params.min_median_read_depth}"
    export SECAT_CHANGEPOINT_METHOD="${params.changepoint_penalty_method}"
    export SECAT_CHANGEPOINT_MULTIPLIER="${params.changepoint_penalty_multiplier}"
    export SECAT_NULL_MODEL_P="${params.null_model_p_threshold}"
    export SECAT_NULL_MODEL_MIN_CONSEC="${params.null_model_min_consecutive}"
    export SECAT_NULL_MODEL_MIN_TRIM_BP="${params.null_model_min_trim_bp}"
    export SECAT_DISTANCE_CUTOFF="${params.distance_cutoff_threshold}"
    export SECAT_DISTANCE_CUTOFF_MIN_BP="${params.distance_cutoff_min_trim_bp}"
    export SECAT_OUTDIR="."
    mkdir -p output/intermediate output/real_data_results/${study_name}
    cp ${study_coords}   output/intermediate/study_alignment_coords.csv 2>/dev/null || true
    cp ${consensus_info} output/intermediate/consensusregioninfo.csv    2>/dev/null || true
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/06_analyse_real.R
    mv output/real_data_results/${study_name}/${study_name}_results.rds ./${study_name}_results.rds
    """
}
