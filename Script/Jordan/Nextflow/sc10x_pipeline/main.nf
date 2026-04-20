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

include { BCL_TO_MULTI        } from './subworkflows/bcl_to_multi'
include { FASTQ_TO_MULTI      } from './subworkflows/fastq_to_multi'

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
    ch_versions = Channel.empty()

    // -----------------------------------------------------------------------
    // STEP 1: Preprocess the sample sheet to standardize it for downstream tools
    // -----------------------------------------------------------------------
    log.info "Step 1: Preprocessing/Standardizing sample sheet."
    
    // Run PREPROCESSING module 
    PREPROCESSING(
        raw_sample_sheet: file(params.sample_sheet, checkIfExists: true), // Access the raw sample sheet path from params
    )

    // Capture outputs from PREPROCESSING
    ch_standardized_sheet = PREPROCESSING.out.standardized_sheet // Capture channel for standardized sample sheet
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions) // Capture channel for versions information (mix for cumulation across steps)

    // -----------------------------------------------------------------------
    // STEP 2: Process BCL files to FASTQ according standardized sample sheet
    // -----------------------------------------------------------------------
    log.info "Step 2: Processing BCL files to FASTQ."
    
    // Run CELLRANGER_MKFASTQ module 
    CELLRANGER_MKFASTQ(
        run_id: params.run_id, // Access the raw sample sheet path from params
        bcl_dir: file(params.input_dir, checkIfExists: true), // BCL directory from params
        sample_sheet: ch_standardized_sheet // Use the standardized sample sheet from the previous step
    )

    // Capture outputs from CELLRANGER_MKFASTQ
    ch_fastq_files = CELLRANGER_MKFASTQ.out.fastq_files // Capture channel for FASTQ files
    ch_versions = ch_versions.mix(CELLRANGER_MKFASTQ.out.versions) // Capture channel for versions information (mix for cumulation across steps)














    // -----------------------------------------------------------------------
    // CASE 1: BCL INPUT -> mkfastq -> multi -> multiqc
    // -----------------------------------------------------------------------
    if (params.input_type == 'bcl') {

        log.info "Mode: BCL detected. Running mkfastq then cellranger multi."

        // Channel: BCL directory + sample sheet
        ch_bcl_input = Channel.of([
            file(params.input_dir, checkIfExists: true),
            file(params.sample_sheet, checkIfExists: true)
        ])

        // Sub-workflow BCL -> MULTI
        BCL_TO_MULTI(ch_bcl_input)

        ch_versions = ch_versions.mix(BCL_TO_MULTI.out.versions)

        // Collect metrics for MultiQC
        ch_qc_inputs = BCL_TO_MULTI.out.metrics
            .collect()
            .map { files -> files }

        MULTIQC(ch_qc_inputs)

    // -----------------------------------------------------------------------
    // CASE 2: FASTQ INPUT -> count -> multiqc
    // -----------------------------------------------------------------------
    } else if (params.input_type == 'fastq') {

        log.info "Mode: FASTQ detected. Auto-detecting samples."

        // Auto-detect FASTQ directories organized by sample
        // Convention: <input_dir>/<sample_id>/*_R1_*.fastq.gz
        ch_fastq_dirs = Channel
            .fromPath("${params.input_dir}/*", type: 'dir')
            .map { dir ->
                def sample_id = dir.name
                def fastq_files = file("${dir}/*.fastq.gz")
                if (fastq_files.size() == 0) {
                    log.warn "WARNING: No FASTQ files found in ${dir} - sample skipped."
                    return null
                }
                return tuple(sample_id, dir)
            }
            .filter { it != null }

        // Ensure at least one sample was found
        ch_fastq_dirs
            .ifEmpty {
                error "ERROR: No FASTQ sample directory found in ${params.input_dir}. Check the directory structure."
            }

        // Sub-workflow FASTQ -> COUNT
        FASTQ_TO_COUNT(ch_fastq_dirs)

        ch_versions = ch_versions.mix(FASTQ_TO_COUNT.out.versions)

        // Collect metrics for MultiQC
        ch_qc_inputs = FASTQ_TO_COUNT.out.metrics
            .collect()
            .map { files -> files }

        MULTIQC(ch_qc_inputs)
    }

    // Emit versions for reproducibility
    ch_versions
        .unique()
        .collectFile(name: "${params.output_dir}/pipeline_versions.txt", newLine: true)

    log.info "✔ Pipeline completed. Results in: ${params.output_dir}"
}

// ========================================================================================
// GLOBAL ERROR HANDLING
// ========================================================================================

workflow.onError {
    log.error """
    ╔══════════════════════════════════════════════════════╗
    ║      PIPELINE FAILED - Check logs for details        ║
    ╚══════════════════════════════════════════════════════╝
    Error  : ${workflow.errorMessage}
    Report : ${workflow.launchDir}/pipeline_report.html
    """.stripIndent()
}

workflow.onComplete {
    if (workflow.success) {
        log.info """
    ╔══════════════════════════════════════════════════════╗
    ║            PIPELINE COMPLETED SUCCESSFULLY           ║
    ╚══════════════════════════════════════════════════════╝
    Total duration : ${workflow.duration}
    Results        : ${params.output_dir}
    HTML report    : ${workflow.launchDir}/pipeline_report.html
    Timeline       : ${workflow.launchDir}/pipeline_timeline.html
    DAG            : ${workflow.launchDir}/pipeline_dag.html
        """.stripIndent()
    }
}







// Lunch the workflow : command to run the workflow, to be executed in the terminal
// nextflow run main.nf -params-file ./params.json -profile slurm