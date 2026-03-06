#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SpaceRanger {
    tag "Processing sample: ${sample_id}"
    publishDir "${params.spaceranger_output}", mode: 'copy'

    container 'cumulusprod/spaceranger:4.0.1'

    input:
    val(sample_id)
    path(transcriptome)

    output:
    tuple val(sample_id), path("${sample_id}_outputs/outs/"), emit: spaceranger_output
    path("${sample_id}_outputs/outs/web_summary.html"), emit: web_summary

    script:
    """
    spaceranger count \\
        --id=${sample_id}_outputs \\
        --transcriptome=${transcriptome} \\
        --probe-set=${params.spaceranger_probe_set} \\
        --libraries=${params.spaceranger_samples} \\
        --feature-ref=${params.spaceranger_feature_reference} \\
        --cytaimage=${params.spaceranger_cytassist_image} \\
        --loupe-alignment=${params.spaceranger_loupe_alignment} \\
        --create-bam=false \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """
}