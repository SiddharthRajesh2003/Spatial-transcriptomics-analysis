#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FastQC {
    tag "Running FastQC on sample: ${sample_id}"
    publishDir "${params.fastqc_dir}/${sample_id}", mode: 'copy'

    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(sample_id), path(fastq_dir)

    output:
    tuple val(sample_id), path("*.zip"), emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc ${fastq_dir}/*.fastq.gz \\
        --threads ${task.cpus} \\
        --outdir .
    """
}