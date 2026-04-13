process GENERATE_VERDICTS {
    tag "generating verdicts"
    label 'mem_4g'
    publishDir "${params.outdir}/aggregated_data", mode: 'copy'

    input:
    path master_verdict_table

    output:
    path "verdict_data_all_levels.csv", emit: verdict_data
    path "final_trim_verdicts.csv",     emit: final_verdicts, optional: true

    script:
    """
    mkdir -p output/aggregated_data
    cp ${master_verdict_table} output/aggregated_data/master_verdict_table.csv
    export SECAT_OUTDIR="."
    export SECAT_PROJECTDIR="${projectDir}"
    Rscript ${projectDir}/R/10_gen_verdicts.R
    cp output/aggregated_data/verdict_data_all_levels.csv ./verdict_data_all_levels.csv
    cp output/aggregated_data/final_trim_verdicts.csv     ./final_trim_verdicts.csv 2>/dev/null || true
    """
}
