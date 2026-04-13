process VALIDATE {
    tag "validation (multi-tier)"
    label 'mem_48g_cpu16'

    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path combined_fasta
    path combined_otu
    path combined_tax
    path combined_meta
    path asv_mapping
    path aggregated_dir

    output:
    path "outputs/**",          emit: validation_outputs
    path "logs/validation.log", emit: validation_log, optional: true

    script:
    """
    mkdir -p logs outputs

    # Reconstruct directory layout expected by validation_taxon_v2.R
    mkdir -p unaligned_cleaned aligned_trimmed output/intermediate output/standardized_datasets

    # --- unaligned_cleaned: pre-trim data from meta_analysis outputs ---
    # The combined outputs from MERGE_DATASETS are the post-trim dataset.
    # Link them as aligned_trimmed; the script compares against unaligned via
    # the asv_mapping_final and trim_summary files.
    cp ${combined_otu}  aligned_trimmed/feature_table.tsv
    cp ${combined_tax}  aligned_trimmed/taxonomy.tsv
    cp ${combined_meta} aligned_trimmed/metadata.tsv       2>/dev/null || true
    cp ${combined_fasta} aligned_trimmed/sequences.fasta   2>/dev/null || true

    # Pre-trim inputs — sourced from aggregated_dir intermediate files
    # (original feature tables are not re-staged through Nextflow;
    #  validation_taxon_v2.R handles graceful missing-file behaviour per tier)
    cp -r ${aggregated_dir}/* output/ 2>/dev/null || true
    cp ${asv_mapping} asv_mapping_final.tsv 2>/dev/null || true

    export SECAT_OUTDIR="."

    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/validation_taxon_v2.R \
        0 \
        FALSE \
        2>&1 | tee logs/validation.log
    """
}
