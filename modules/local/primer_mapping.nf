process PRIMER_MAPPING {
    tag "primer coordinate mapping"
    label 'mem_64g'
    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path reference_db

    output:
    path "primer_coords_phase1_output.csv", emit: primer_coords

    script:
    """
    export SECAT_REFERENCE_DB="${reference_db}"
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_MAX_PRIMER_MISMATCH="${params.max_primer_mismatch}"
    export SECAT_OUTDIR="."
    mkdir -p output/intermediate
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/01_primer_mapping.R
    cp output/intermediate/primer_coords_phase1_output.csv ./primer_coords_phase1_output.csv
    """
}

process GENERATE_PRIMER_DBS {
    tag "${primer_name}"
    label 'mem_16g'
    publishDir "${params.outdir}/primer_databases", mode: 'copy', pattern: "db_*.fasta"

    input:
    tuple val(primer_name), path(reference_db)

    output:
    path "db_${primer_name}.fasta", emit: primer_db

    script:
    """
    export SECAT_REFERENCE_DB="${reference_db}"
    export SECAT_OUTDIR="."
    mkdir -p output/intermediate output/primer_databases
    cp ${params.outdir}/intermediate/primer_coords_phase1_output.csv \
       output/intermediate/primer_coords_phase1_output.csv
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/00_generate_primers.R "${primer_name}"
    cp output/primer_databases/db_${primer_name}.fasta ./db_${primer_name}.fasta
    """
}
