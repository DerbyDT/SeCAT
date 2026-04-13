process GENERATE_INDEX {
    tag "generating HTML index"
    label 'mem_8g'
    publishDir "${params.outdir}/reports", mode: 'copy', pattern: "*.html"

    input:
    path aggregated_dir
    path verdict_data

    output:
    path "*.html", emit: html_index, optional: true

    script:
    """
    mkdir -p output/aggregated_data output/reports
    cp -r ${aggregated_dir}/* output/aggregated_data/ 2>/dev/null || true
    cp ${verdict_data} output/aggregated_data/verdict_data_all_levels.csv
    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/09_gen_index.R
    cp output/reports/*.html ./ 2>/dev/null || true
    """
}
