process VALIDATE {
    tag "validation (multi-tier)"
    label 'mem_30g'
    publishDir "${params.outdir}/validation", mode: 'copy'

input:
    path post_otu,    stageAs: 'post_feature_table.tsv'
    path post_tax,    stageAs: 'post_taxonomy.tsv'
    path post_meta,   stageAs: 'post_metadata.tsv'
    path post_fasta,  stageAs: 'post_sequences.fasta'
    path pre_otu,     stageAs: 'pre_feature_table.tsv'
    path pre_tax,     stageAs: 'pre_taxonomy.tsv'
    path pre_meta,    stageAs: 'pre_metadata.tsv'
    path asv_mapping
    path aggregated_dir

    output:
    path "outputs/**",          emit: validation_outputs
    path "logs/validation.log", emit: validation_log, optional: true

    script:
    """
    mkdir -p logs outputs pre_consensus post_consensus \
             output/intermediate output/standardized_datasets

    # Stage post-consensus (trimmed MetaASV dataset)
    cp post_feature_table.tsv  post_consensus/feature_table.tsv
    cp post_taxonomy.tsv       post_consensus/taxonomy.tsv
    cp post_metadata.tsv       post_consensus/metadata.tsv   2>/dev/null || true
    cp post_sequences.fasta    post_consensus/sequences.fasta 2>/dev/null || true

    # Stage pre-consensus (original untrimmed dataset)
    cp pre_feature_table.tsv   pre_consensus/feature_table.tsv
    cp pre_taxonomy.tsv        pre_consensus/taxonomy.tsv
    cp pre_metadata.tsv        pre_consensus/metadata.tsv    2>/dev/null || true

    # Stage SeCAT intermediate outputs for coordinate verification (Tier 0C)
    cp -r ${aggregated_dir}/* output/ 2>/dev/null || true
    cp ${asv_mapping} asv_mapping_final.tsv        2>/dev/null || true

    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    export SECAT_VALIDATION_LEVELS="${params.validation_levels ?: "ASV,Genus,Family"}"

    Rscript ${projectDir}/R/validation_taxon_v2.R \
        0 \
        FALSE \
        2>&1 | tee logs/validation.log
    """
}
