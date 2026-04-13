process PREPARE_SIMS {
    tag "preparing simulation tasks"
    label 'mem_4g'
    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path study_coords
    val reference_db
    path clean_manifest

    output:
    path "simulation_tasks.csv",              emit: sim_tasks_csv
    path "simulation_reference_subset.fasta", emit: sim_reference_subset
    path "consensusregioninfo.csv",           emit: consensus_info
    path "study_alignment_coords.csv",        emit: study_align_coords, optional: true

    script:
    """
    export SECAT_REFERENCE_DB="${reference_db}"
    export SECAT_MANIFEST="\$(realpath ${clean_manifest})"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_NUM_SIMULATIONS="${params.num_simulations}"
    export SECAT_SIM_MAX_SILVA_SUBSET="${params.sim_max_silva_subset}"
    export SECAT_SIM_USE_PREBUILT="${params.sim_use_prebuilt_subset}"
    export SECAT_TRIM_STEP_MODE="${params.trim_step_mode}"
    export SECAT_DEFAULT_MAX_TRIM_STEPS="${params.default_max_trim_steps}"
    export SECAT_CONSENSUS_BUFFER_STEPS="${params.consensus_buffer_steps}"
    export SECAT_CONSENSUS_OPT_THRESHOLD="${params.consensus_optimization_threshold}"
    export SECAT_MIN_CONSENSUS_STUDIES="${params.min_consensus_studies}"
    export SECAT_OUTDIR="."
    mkdir -p intermediate
    cp ${study_coords} ./intermediate/study_alignment_coords.csv
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/04_prepare_sims.R
    cp ./intermediate/simulation_tasks.csv              ./simulation_tasks.csv
    cp ./intermediate/simulation_reference_subset.fasta ./simulation_reference_subset.fasta
    cp ./intermediate/consensusregioninfo.csv           ./consensusregioninfo.csv
    cp ./intermediate/study_alignment_coords.csv        ./study_alignment_coords.csv 2>/dev/null || true
    """
}
