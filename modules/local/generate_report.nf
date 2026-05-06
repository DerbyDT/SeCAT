process GENERATE_REPORT {
    tag "generating PDF report"
    label 'mem_64g_cpu4'
    publishDir "${params.outdir}/reports",     mode: 'copy', pattern: "*.pdf"
    publishDir "${params.outdir}/final_plots", mode: 'copy', pattern: "final_plots/**"

    input:
    path aggregated_dir
    path verdict_data
    path rds_files
    path study_coords
    path consensus_info

    output:
    path "*.pdf",          emit: pdf_reports, optional: true
    path "final_plots/**", emit: plot_files,  optional: true

    script:
    """
    mkdir -p output/aggregated_data output/final_plots output/reports \
             output/real_data_results output/intermediate

    cp -r ${aggregated_dir}/* output/aggregated_data/ 2>/dev/null || true
    cp ${verdict_data} output/aggregated_data/verdict_data_all_levels.csv

    for f in ${rds_files}; do
        cp \$f output/real_data_results/
    done

    cp ${study_coords} output/intermediate/study_alignment_coords.csv
    cp ${consensus_info} output/intermediate/consensusregioninfo.csv

    export SECAT_MANIFEST="${params.manifest}"
    export SECAT_ANALYSIS_MODE="${params.analysis_mode}"
    export SECAT_PLOT_DPI="${params.plot_dpi}"
    export SECAT_FORCE_REGENERATE="${params.force_regenerate}"
    export SECAT_OUTDIR="./output"
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/08_gen_report.R

    cp output/reports/*.pdf ./ 2>/dev/null || true
    cp -r output/final_plots ./ 2>/dev/null || true
    """
}
