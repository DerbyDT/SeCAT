process MERGE_DATASETS {
    tag "merging datasets"
    label 'mem_30g'
    publishDir "${params.outdir}/meta_analysis", mode: 'copy'

    input:
    path standardized_fastas
    path trim_summary
    path manifest_rows

    output:
    path "combined_sequences.fasta",   emit: combined_fasta
    path "combined_feature_table.tsv", emit: combined_otu
    path "combined_taxonomy.tsv",      emit: combined_tax
    path "combined_metadata.tsv",      emit: combined_meta,  optional: true
    path "asv_mapping_final.tsv",      emit: asv_mapping,    optional: true

    script:
    """
    mkdir -p output/standardized_datasets output/aggregated_data output/meta_analysis output/intermediate
    for f in ${standardized_fastas}; do cp \$f output/standardized_datasets/; done
    cp ${trim_summary} output/aggregated_data/trim_summary.csv
    ln -s ${params.outdir}/intermediate/aligned_fastas output/intermediate/aligned_fastas 2>/dev/null || true
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_VSEARCH_PATH="vsearch"  # binary is on PATH inside container
    export SECAT_MERGE_METHOD="${params.merge_method}"
    export SECAT_HARMONIZE_METADATA="${params.harmonize_metadata}"
    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/13_merge_datasets.R
    cp output/meta_analysis/combined_sequences.fasta      ./ 2>/dev/null || true
    cp output/meta_analysis/combined_feature_table.tsv    ./ 2>/dev/null || true
    cp output/meta_analysis/combined_taxonomy.tsv         ./ 2>/dev/null || true
    cp output/meta_analysis/combined_metadata.tsv         ./ 2>/dev/null || true
    cp output/meta_analysis/asv_mapping_final.tsv         ./ 2>/dev/null || true
    """
}
