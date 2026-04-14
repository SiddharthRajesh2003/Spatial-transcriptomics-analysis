#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CellBender {
    tag "Running CellBender on ${sample_id}"
    publishDir "${params.cellbender_output}", mode: 'copy'

    container 'cumulusprod/cellbender:0.3.0'

    input:
    tuple val(sample_id), path(raw_h5)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_cellbender.h5"), emit: cellbender_h5
    path("${sample_id}/"), emit: cellbender_dir

    script:
    """
    mkdir -p ${sample_id}
    cellbender remove-background \\
        --input ${raw_h5} \\
        --output ${sample_id}/${sample_id}_cellbender.h5 \\
        --expected-cells ${params.cellbender_expected_cells} \\
        --total-droplets-included ${params.cellbender_total_droplets} \\
        --epochs ${params.cellbender_epochs} \\
        --cuda
    """
}
