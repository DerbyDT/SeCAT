process AGGREGATE {
    tag "aggregating all results"
    label 'mem_64g'
    publishDir "${params.outdir}/aggregated_data", mode: 'copy', pattern: "aggregated_data/*"

    input:
    path real_results
    path sim_results
    path study_coords
    path consensus_info

    output:
    path "aggregated_data/master_verdict_table.csv",           emit: master_verdict_table
    path "aggregated_data/verdict_data_all_levels.csv",        emit: verdict_data,     optional: true
    path "aggregated_data/simulation_baseline_statistics.csv", emit: baseline_stats,   optional: true
    path "aggregated_data/simulation_retention_curves.csv",    emit: retention_curves, optional: true
    path "aggregated_data",                                    emit: aggregated_dir

    script:
    """
    mkdir -p output/intermediate output/real_data_results output/simulation_results aggregated_data
    cp ${study_coords}   output/intermediate/study_alignment_coords.csv 2>/dev/null || true
    cp ${consensus_info} output/intermediate/consensusregioninfo.csv    2>/dev/null || true
    for f in *_results.rds; do
        [[ "\$f" == *__seed_* ]] && continue
        study=\$(basename \$f _results.rds)
        mkdir -p output/real_data_results/\${study}
        cp \$f output/real_data_results/\${study}/\${study}_results.rds
    done
    for f in *__seed_*__results.rds; do
        fname=\$(basename \$f)
        task_id=\$(echo \$fname | sed 's/__seed_.*//')
        seed=\$(echo \$fname | sed 's/.*__seed_//' | sed 's/__results.rds//')
        mkdir -p output/simulation_results/\${task_id}/seed_\${seed}
        cp \$f output/simulation_results/\${task_id}/seed_\${seed}/results.rds
    done
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_CHANGEPOINT_METHOD="${params.changepoint_penalty_method}"
    export SECAT_CHANGEPOINT_MULTIPLIER="${params.changepoint_penalty_multiplier}"
    export SECAT_NULL_MODEL_P="${params.null_model_p_threshold}"
    export SECAT_NULL_MODEL_MIN_CONSEC="${params.null_model_min_consecutive}"
    export SECAT_DISTANCE_CUTOFF="${params.distance_cutoff_threshold}"
    export SECAT_CONSENSUS_OPT_THRESHOLD="${params.consensus_optimization_threshold}"
    export SECAT_MIN_CONSENSUS_STUDIES="${params.min_consensus_studies}"
    export SECAT_OUTDIR="./output"
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/07_aggregate.R
    cp output/aggregated_data/* aggregated_data/ 2>/dev/null || true
    """
}
