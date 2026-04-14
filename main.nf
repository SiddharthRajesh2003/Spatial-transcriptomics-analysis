#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Import modules
include { FastQC } from './modules/fastqc/'
include { CellRanger } from './modules/cellranger/'
include { CellBender } from './modules/cellbender/'
include { SpaceRanger } from './modules/spaceranger/'
include { MultiQC } from './modules/multiqc/'

def helpMessage() {
    log.info"""
    =======================================================================
                Spatial and Single Cell Processing Pipeline
    =======================================================================

    Nextflow pipeline for processing spatial and single cell RNA-seq data.
    Leverages CellRanger, SpaceRanger, FastQC, and MultiQC for comprehensive analysis.

    Usage:
        nextflow run main.nf [options]
    
    Required Options:
        --cellranger_samples    Path to samplesheet CSV (columns: sample_id, fastq_dir)
        --spaceranger_samples   Path to libraries CSV for SpaceRanger (columns: fastqs, sample, library_type)

    Optional Parameters:
        --results               Output directory [default: \${base}/results]
        --fastqc_dir            FastQC output directory [default: \${results}/fastqc]
        --cellranger_output     Cell Ranger output directory [default: \${results}/cellranger]
        --spaceranger_output    Space Ranger output directory [default: \${results}/spaceranger]
        --transcriptome         Path to reference transcriptome [default: set in nextflow.config]

    Skip Options:
        --skip_fastqc           Skip FastQC step [default: false]
        --skip_cellranger       Skip Cell Ranger and use existing outputs [default: false]
        --skip_spaceranger      Skip Space Ranger and use existing outputs [default: false]
        --fallback_to_cellranger    Run Cell Ranger if existing outputs not found [default: true]
        --fallback_to_spaceranger   Run Space Ranger if existing outputs not found [default: true]

    Example:
        nextflow run main.nf -profile slurm \\
            --cellranger_samples samplesheets/cellranger_samples.csv \\
            --spaceranger_samples samplesheets/spaceranger_samples.csv
    """
}

// Validate required parameters
def validateParams() {
    def errors = []
    
    if (!params.cellranger_samples) {
        errors << "Missing required parameter: --cellranger_samples"
    }
    if (!params.spaceranger_samples) {
        errors << "Missing required parameter: --spaceranger_samples"
    }
    if (!params.transcriptome) {
        errors << "Missing required parameter: --transcriptome"
    }    
    if (errors) {
        log.error "Validation errors:"
        errors.each { error -> log.error "  - ${error}" }
        log.error "\nRun 'nextflow run main.nf --help' for usage information"
        exit 1
    }
}

// Check if Cell Ranger should actually be skipped
def actuallyskipCellRanger() {
    if (!params.skip_cellranger) {
        return false
    }

    def cellrangerDir = file(params.cellranger_output)
    if (!cellrangerDir.exists()) {
        if (params.fallback_to_cellranger) {
            log.warn "Cell Ranger directory ${params.cellranger_output} does not exist! Will run cellranger!"
            return false
        } else {
            error "No Cell Ranger outputs found in directory ${params.cellranger_output} and fallback_to_cellranger disabled!"
        }
    }
    return true
}

def actuallyskipSpaceRanger(){
    if (!params.skip_spaceranger) {
        return false
    }

    def spacerangerDir = file(params.spaceranger_output)
    if (!spacerangerDir.exists()){
        if (params.fallback_to_spaceranger) {
            log.warn "Space Ranger directory ${params.spaceranger_output} does not exist! Will run spaceranger!"
            return false
        } else {
            error "No Space Ranger outputs found in directory ${params.spaceranger_output} and fallback_to_spaceranger disabled!"
        }
    }
    return true
}

// Main workflow
workflow {
    if (params.help){
        helpMessage()
        exit 0
    }

    // Validate parameters
    validateParams()

    // Read CellRanger and SpaceRanger sample sheets
    spaceranger_qc_samples_ch = channel
                                    .fromPath(params.spaceranger_samples, checkIfExists: true)
                                    .splitCsv(header: true)
                                    .map { row -> tuple(row.sample, file(row.fastqs)) }

    cellranger_qc_samples_ch = channel
                                    .fromPath(params.cellranger_samples, checkIfExists: true)
                                    .splitCsv(header: true)
                                    .map { row -> tuple(row.sample_id, file(row.fastq_dir)) }
    
    cellranger_samples_ch = channel
                                .fromPath(params.cellranger_samples, checkIfExists: true)
                                .splitCsv(header: true)
                                .map{
                                    row ->
                                    tuple(row.sample_id, file("${row.fastq_dir}"))
                                }
    
    spaceranger_sample_id_ch = channel.value(params.spaceranger_id)


    if (!params.skip_fastqc){
        log.info "Running FastQC on all samples..."
        FastQC(cellranger_qc_samples_ch.mix(spaceranger_qc_samples_ch))
        qc_reports = FastQC.out.zip.map { _id, zip -> zip }.collect()
    } else {
        log.info "Skipping FastQC step as per user request..."
        qc_reports = channel.value([])
    }
    
    // Determine whether to run Cell Ranger or use existing outputs
    def skipCellRanger = actuallyskipCellRanger()

    if (skipCellRanger) {
        log.info "Using existing Cell Ranger outputs from ${params.cellranger_output}"

        // Create channel from existing Cell Ranger outputs
        cellranger_outputs_ch = channel
            .fromPath("${params.cellranger_output}/*/outs/filtered_feature_bc_matrix", type: 'dir', checkIfExists: true)
            .map { matrix_dir ->
                // Extract sample_id from directory structure
                // Assumes structure: cellranger_output/sample_id_output/outs/filtered_feature_bc_matrix
                def sample_id = matrix_dir.parent.parent.name.replaceAll("_output", "")
                tuple(sample_id, matrix_dir)
            }

        // Create channel for existing sample output directories (preserves unique sample names)
        cellranger_reports = channel
            .fromPath("${params.cellranger_output}/*_output", type: 'dir', checkIfExists: true)
            .toList()

    } else {
        log.info "Running Cell Ranger on samples..."

        // Run Cell Ranger
        CellRanger(
            cellranger_samples_ch,
            params.transcriptome
        )

        // Extract filtered_feature_bc_matrix for downstream analysis
        cellranger_outputs_ch = CellRanger.out.cellranger_out.map {
            sample_id, outs_dir ->
                tuple(sample_id, file("${outs_dir}/filtered_feature_bc_matrix"))
        }

        // Collect sample output directories for MultiQC (preserves unique sample names)
        cellranger_reports = CellRanger.out.cellranger_out
            .map { _sample_id, outs_dir -> outs_dir.parent }
            .toList()
        
        CellBender(
            CellRanger.out.raw_h5    
        )
    }

    def skipSpaceRanger = actuallyskipSpaceRanger()

    if (skipSpaceRanger) {
        log.info "Using existing Space Ranger outputs from ${params.spaceranger_output}"

        // Create channel from existing Space Ranger outputs
        spaceranger_outputs_ch = channel
            .fromPath("${params.spaceranger_output}/*/outs/filtered_feature_bc_matrix", type: 'dir', checkIfExists: true)
            .map { matrix_dir ->
                // Extract sample_id from directory structure
                def sample_id = matrix_dir.parent.parent.name.replaceAll("_outputs", "")
                tuple(sample_id, matrix_dir)
            }
        
        // Create channel for existing sample output directories (preserves unique sample names)
        spaceranger_reports = channel
            .fromPath("${params.spaceranger_output}/*_outputs", type: 'dir', checkIfExists: true)
            .toList()
    } else {
        log.info "Running Space Ranger on samples..."

        // Run Space Ranger
        SpaceRanger(spaceranger_sample_id_ch, params.transcriptome)

        // Extract filtered_feature_bc_matrix for downstream analysis
        spaceranger_outputs_ch = SpaceRanger.out.spaceranger_output.map {
            sample_id, outs_dir ->
                tuple(sample_id, file("${outs_dir}/filtered_feature_bc_matrix"))
        }

        spaceranger_reports = SpaceRanger.out.spaceranger_output
            .map { _sample_id, outs_dir -> outs_dir.parent }
            .toList()
    }

    if (!params.skip_fastqc || !params.skip_cellranger || !params.skip_spaceranger) {
        log.info "Running MultiQC to aggregate QC reports..."
        MultiQC(
            qc_reports,
            cellranger_reports,
            spaceranger_reports
        )
    }
}