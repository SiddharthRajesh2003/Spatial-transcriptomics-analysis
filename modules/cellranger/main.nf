#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CellRanger {
    tag "Running CellRanger count for sample: ${sample_id}"
    publishDir "${params.cellranger_output}", mode: 'copy'

    container 'cumulusprod/cellranger:10.0.0'

    input:
    tuple val(sample_id), path(fastq_dir)
    path transcriptome
    
    output:
    tuple val(sample_id), path("${sample_id}_outputs/outs/"), emit: cellranger_out
    tuple val(sample_id), path ()"${sample_id}_outputs/outs/raw_feature_bc_matrix.h5"), emit: raw_h5
    path("${sample_id}_outputs/outs/web_summary.html"), emit: web_summary

    script:
    """
    cellranger count \\
        --id=${sample_id}_outputs \\
        --transcriptome=${transcriptome} \\
        --fastqs=${fastq_dir} \\
        --sample=${sample_id} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --create-bam=false
    """
}