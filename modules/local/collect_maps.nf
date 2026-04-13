process COLLECT_MAPS {
    tag "collecting ${mapping_parts.size()} studies"
    label 'mem_4g'
    publishDir "${params.outdir}/intermediate", mode: 'copy'

    input:
    path mapping_parts

    output:
    path "study_alignment_coords.csv",  emit: study_coords
    path "study_mapping_summary.csv",   emit: collected_maps, optional: true

    script:
    """
    # Collect all mapping part files passed in by Nextflow
    PARTS=(\$(ls mapping_part_*.csv 2>/dev/null | sort -t_ -k3 -n))

    if [ \${#PARTS[@]} -eq 0 ]; then
        echo "ERROR: No mapping_part_*.csv files found in work directory"
        exit 1
    fi

    echo "Collecting \${#PARTS[@]} mapping parts..."

    # Write header from first file
    head -1 "\${PARTS[0]}" > study_alignment_coords.csv

    # Append data rows from all parts in order
    for f in "\${PARTS[@]}"; do
        tail -n +2 "\$f" >> study_alignment_coords.csv
    done

    NROWS=\$(tail -n +2 study_alignment_coords.csv | wc -l)
    echo "Aggregated \${NROWS} study coordinate rows"

    # Summary file for downstream reference
    echo "studies_collected,\${#PARTS[@]}" > study_mapping_summary.csv
    echo "coordinate_rows,\${NROWS}" >> study_mapping_summary.csv
    """
}
