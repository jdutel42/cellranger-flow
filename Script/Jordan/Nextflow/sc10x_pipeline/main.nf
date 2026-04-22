#!/usr/bin/env nextflow

/*
========================================================================================
    Cellranger_Nextflow : Single-Cell 10x Genomics Pipeline (DSL2)
========================================================================================
    Author    : Jordan Dutel - GENIM Team, CRCT
    Version   : 1.0.0
    Usage     : nextflow run main.nf -profile docker --input_type fastq --input_dir /data/fastq
----------------------------------------------------------------------------------------
*/

// Activate DSL2 language, which is required for modular workflows and subworkflows
nextflow.enable.dsl = 2

// ========================================================================================
// IMPORTS MODULS
// ========================================================================================

include { PREPROCESSING       } from './modules/pre_processing'
include { CELLRANGER_MKFASTQ  } from './modules/cellranger_mkfastq'
include { CELLRANGER_MULTI    } from './modules/cellranger_multi'
include { MULTIQC             } from './modules/multiqc'

// ========================================================================================
// IMPORTS SUB-WORKFLOWS
// ========================================================================================



// ========================================================================================
// DISPLAY INITIAL LOGGING
// ========================================================================================

def printHeader() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║         Cellranger_Nextflow — Single-Cell 10x Genomics DSL2      ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  input_type    : ${params.input_type}                            ║
    ║  input_dir     : ${params.input_dir}                             ║
    ║  output_dir    : ${params.output_dir}                            ║
    ║  genome        : ${params.genome_reference}                      ║
    ║  sample_sheet  : ${params.sample_sheet}                          ║
    ║  run_id        : ${params.run_id}                                ║
    ║  localcores    : ${params.localcores}                            ║
    ║  localmemory   : ${params.localmemory} GB                        ║
    ╚══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

// ========================================================================================
// REQUIRED PARAMETERS VALIDATION
// ========================================================================================

def validateParams() {

    // Input type validation
    if (!params.input_type) {
        error "ERROR: The --input_type parameter is required ('bcl' or 'fastq')."
    }
    if (!['bcl', 'fastq'].contains(params.input_type)) {
        error "ERROR: --input_type must be 'bcl' or 'fastq'. Received value: '${params.input_type}'."
    }

    // Input directory validation
    if (!params.input_dir) {
        error "ERROR: The --input_dir parameter is required."
    }
    def input_path = file(params.input_dir)
    if (!input_path.exists()) {
        error "ERROR: Input directory does not exist: ${params.input_dir}"
    }
    if (!input_path.isDirectory()) {
        error "ERROR: --input_dir must be a directory, not a file: ${params.input_dir}"
    }

    // Genome reference validation
    if (!params.genome_reference) {
        error "ERROR: The --genome_reference parameter is required."
    }
    def genome_path = file(params.genome_reference)
    if (!genome_path.exists()) {
        error "ERROR: Genome reference directory does not exist: ${params.genome_reference}"
    }

    // Sample sheet validation for BCL mode
    if (params.input_type == 'bcl') {
        if (!params.sample_sheet) {
            error "ERROR: --sample_sheet is required in BCL mode."
        }
        def ss_path = file(params.sample_sheet)
        if (!ss_path.exists()) {
            error "ERROR: Sample sheet does not exist: ${params.sample_sheet}"
        }
    }

    // Output directory validation
    if (!params.output_dir) {
        error "ERROR: The --output_dir parameter is required."
    }

    log.info "✔ Parameter validation passed."
}

// ========================================================================================
// MAIN WORKFLOW
// ========================================================================================

workflow {

    // Print header and validate parameters at the start of the workflow
    printHeader()

    // Validate parameters before starting any processing
    validateParams()

    // Channel for version files (for MultiQC and traceability)
    // Initialization is needed to ensure accumulation versions from multiple steps
    ch_versions = channel.empty()

    // -----------------------------------------------------------------------
    // STEP 1: Preprocess the sample sheet to standardize it for downstream tools
    // -----------------------------------------------------------------------
    log.info "Step 1: Preprocessing/Standardizing sample sheet."
    
    // Run PREPROCESSING module 
    PREPROCESSING(
        raw_sample_sheet_file_path: file(params.raw_sample_sheet_file_path, checkIfExists: true), // Access the raw sample sheet path from params
        run_id: params.run_id, // Pass run_id for logging and output naming
        preprocessing_output_dir: params.preprocessing_output_dir, // Pass output directory for standardized sample sheet
        preprocessing_log_dir: params.preprocessing_log_dir // Pass log directory for preprocessing logs
    )

    // Capture outputs from PREPROCESSING
    ch_preprocessed_sample_sheet = PREPROCESSING.out.preprocessed_sample_sheet // Capture channel for standardized sample sheet
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions) // Capture channel for versions information (mix for cumulation across steps)

    log.info "✔ Sample sheet preprocessing completed. Standardized sample sheet available at: ${params.preprocessing_output_dir}"

    // -----------------------------------------------------------------------
    // STEP 2: Process BCL files to FASTQ according standardized sample sheet
    // -----------------------------------------------------------------------
    log.info "Step 2: Processing BCL files to FASTQ."
    
    // Run CELLRANGER_MKFASTQ module 
    CELLRANGER_MKFASTQ(
        run_id: params.run_id, // Access the raw sample sheet path from params
        bcl_dir: params.bcl_dir, // BCL directory from params
        preprocessed_sample_sheet: ch_preprocessed_sample_sheet, // Use the standardized sample sheet from the previous step
        qc_output_dir: params.qc_output_dir, // Pass output directory for FASTQ files
        qc_log_dir: params.qc_log_dir // Pass log directory for QC logs
    )

    // Capture outputs from CELLRANGER_MKFASTQ
    ch_fastqs = CELLRANGER_MKFASTQ.out.fastqs // Capture channel for generated FASTQ files
    ch_versions = ch_versions.mix(CELLRANGER_MKFASTQ.out.versions) // Capture channel for versions information (mix for cumulation across steps)
    
    log.info "✔ BCL to FASTQ processing completed. FASTQ files available at: ${params.qc_output_dir}"

    // -----------------------------------------------------------------------
    // STEP 3: Perform Alignment with Cellranger Multi
    // -----------------------------------------------------------------------
    log.info "Step 3: Performing alignment with Cellranger Multi."

    // Convert the comma-separated batch IDs string into a list of batch id
    batch_ids_list = params.batch_ids.toString().split(',') // Split the comma-separated batch IDs into a list of batch_id

    // Create a channel with several elements (batch_id), each will be pass to a separate instance of CELLRANGER_MULTI, enabling parallel alignment across all requested batches.
    ch_batch_id = channel
        .from(batch_ids_list) // From the list of batch IDs
        .map { batch_id -> "${params.protocol_prefix}_batch${batch_id}" } // Map each batch_id in batch_ids_list to a full_batch_name
    
    // Run CELLRANGER_MULTI module
    CELLRANGER_MULTI(
        run_id: params.run_id, // Pass run_id for logging and output naming
        batch_id: ch_batch_id, // Pass the channel of batch identifiers
        fastq_folder: ch_fastqs, // Pass the channel of FASTQ directories generated from the previous step
        fastqs: ch_fastqs, // Pass the channel of FASTQ directories generated from the previous step
        ref_gex: file(params.path_ref_gex, checkIfExists: true), // GEX reference from params
        ref_vdj: file(params.path_ref_vdj, checkIfExists: true), // VDJ reference from params
        alignment_output_dir: params.alignment_output_dir, // Pass output directory for Cellranger Multi results
        alignment_log_dir: params.alignment_log_dir // Pass log directory for Cellranger Multi
    )

    // Capture outputs from CELLRANGER_MULTI
    ch_metrics = CELLRANGER_MULTI.out.metrics // Capture channel for metrics summary
    ch_web_summaries = CELLRANGER_MULTI.out.web_summaries // Capture channel for web summaries
    ch_versions = ch_versions.mix(CELLRANGER_MULTI.out.versions) // Capture channel for versions information (mix for cumulation across steps)

    log.info "✔ Cellranger Multi processing completed. Results available at: ${params.alignment_output_dir}"

    // -----------------------------------------------------------------------
    // STEP 4: MultiQC report generation
    // -----------------------------------------------------------------------
    log.info "Step 4: Generating MultiQC report."
    
    // Combine metrics summaries and web summaries into a single channel for MultiQC input
    ch_qc_files = 
        ch_metrics.map { _batch_id, f -> f } // In ch_metrics (that is a tuple(batch_id, metrics_file_path)), extract the file path (f) for each batch_id and create a channel of metrics summary files
        .mix(ch_web_summaries.map { _batch_id, f -> f }) 
        // In ch_web_summaries (that is a tuple(batch_id, web_summary_file_path)), 
        // extract the file path (f) for each batch_id and create a channel of web summary files 
        // ==> then mix it with the channel of metrics summary files to create a single channel of QC files for MultiQC input, containing both metrics summaries and web summaries from all batches

    // Run MULTIQC module
    MULTIQC(
        run_id: params.run_id, // Pass run_id for logging and output naming
        qc_files: ch_qc_files.collect(), // Pass the channel of QC files collected (metrics summaries and web summaries) generated from Cellranger Multi
        multiqc_output_dir: params.multiqc_output_dir, // Pass output directory for MultiQC results
        multiqc_log_dir: params.multiqc_log_dir // Pass log directory for MultiQC logs
    )


    
    
    
    
}    
// ========================================================================================
// GLOBAL ERROR HANDLING
// ========================================================================================







// Lunch the workflow : command to run the workflow, to be executed in the terminal
// nextflow run main.nf -params-file ./params.json -profile slurm -with-conda