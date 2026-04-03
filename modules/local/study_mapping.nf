process STUDY_MAPPING {
    tag "${study_name}"
    label 'mem_64g_cpu4'
    publishDir "${params.outdir}/intermediate/study_mapping_parts", mode: 'copy', pattern: "mapping_part_*.csv"
    publishDir "${params.outdir}/intermediate/asv_coordinates",     mode: 'copy', pattern: "*_coords.csv"
    publishDir "${params.outdir}/intermediate/aligned_fastas",      mode: 'copy', pattern: "*_aligned.fasta"

    input:
    tuple val(task_id), val(study_name), val(manifest_row)
    path reference_db

    output:
    path "mapping_part_${task_id}.csv",  emit: mapping_part
    path "${study_name}_coords.csv",     emit: asv_coords,    optional: true
    path "${study_name}_aligned.fasta",  emit: aligned_fasta, optional: true

    script:
    """
    export SGE_TASK_ID=${task_id}
    export SECAT_REFERENCE_DB="${reference_db}"
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_ALIGNMENT_MODE="${params.reference_alignment_mode}"
    export SECAT_SUBSET_SIZE="${params.reference_subset_size}"
    export SECAT_USE_ALL_ASVS="${params.use_all_asvs_for_mafft}"
    export SECAT_VSEARCH_PATH="vsearch"  # binary is on PATH inside container
    export SECAT_ALIGNMENT_METHOD="${params.study_alignment_method}"
    export SECAT_OUTDIR="${params.outdir}"
    Rscript ${projectDir}/R/02_study_mapping.R
    """
}
