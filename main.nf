#!/usr/bin/env nextflow

/*
========================================================================================
    Cellranger_Nextflow : Single-Cell 10x Genomics Pipeline (DSL2)
========================================================================================
    Author    : Jordan Dutel - GENIM Team, CRCT
    Version   : 1.0.0
    Usage     : nextflow run main.nf -params-file ./params.json -profile slurm -with-conda
----------------------------------------------------------------------------------------
*/

// Activate DSL2 language, which is required for modular workflows and subworkflows
nextflow.enable.dsl = 2

// ========================================================================================
// IMPORTS MODULS
// ========================================================================================

include { PREPROCESSING       } from './modules/preprocessing.nf'
include { CELLRANGER_MKFASTQ  } from './modules/cellranger_mkfastq.nf'
include { CELLRANGER_MULTI    } from './modules/cellranger_multi.nf'
include { MULTIQC             } from './modules/multiqc.nf'
include { MERGE_VERSIONS      } from './modules/merge_versions.nf'

// ========================================================================================
// IMPORTS SUB-WORKFLOWS
// ========================================================================================

// (No sub-workflows in this version, but this section is reserved for future additions if needed)

// ========================================================================================
// DISPLAY INITIAL LOGGING
// ========================================================================================

def printHeader() {
    log.info """
    ╔════════════════════════════════════════════════════════════════════════════╗
    ║              CellRanger_Flow — Single-Cell 10x Genomics DSL2               ║
    ╠════════════════════════════════════════════════════════════════════════════╣
    ║  run_id                   : ${params.run_id}
    ║  raw_sample_sheet         : ${params.raw_sample_sheet_file_path}
    ║  bcl_dir                  : ${params.bcl_dir}
    ║  protocol_prefix          : ${params.protocol_prefix}
    ║  batch_ids                : ${params.batch_ids}
    ║  ref_gex                  : ${params.path_ref_gex}
    ║  ref_vdj                  : ${params.path_ref_vdj}
    ║  output_dir               : ${params.output_dir}
    ║  preprocessing_output_dir : ${params.preprocessing_output_dir}
    ║  qc_output_dir            : ${params.qc_output_dir}
    ║  alignment_output_dir     : ${params.alignment_output_dir}
    ║  multiqc_output_dir       : ${params.multiqc_output_dir}
    ║  log_dir                  : ${params.log_dir}
    ║  cpu/memory_limit         : ${params.cpu_limit} / ${params.memory_limit}GB
    ║  pipeline_version         : ${params.pipeline_version}
    ║  min_nextflow_version     : ${params.min_nextflow_version}
    ╚════════════════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

// ========================================================================================
// REQUIRED PARAMETERS VALIDATION
// ========================================================================================

def validateParams() {

    // Log the start of parameter validation
    log.info "[INFO]: 🔍 Validating input parameters..."

    // Check for minimum Nextflow version
    if (params.min_nextflow_version) {
        def current_version = nextflow.version
        if (current_version < params.min_nextflow_version) {
            error "[ERROR]: ❌ Nextflow version ${params.min_nextflow_version} or higher is required. Current version: ${current_version}"
        }
    }

    // Check for pipeline version (optional, but can be used for traceability)
    if (!params.pipeline_version) {
        log.warn "[WARNING]: ⚠️ --pipeline_version is not specified. It's recommended to provide a version for traceability."
    }

    // Check for run_id (required for logging and output naming)
    if (!params.run_id) {
        error "[ERROR]: ❌ The --run_id parameter is required for logging and output naming. Please provide a unique identifier for this run like Flowcell ID, e.g., 'HCHNTDMX2'."
    }

    // Check for raw sample sheet file path
    if (!params.raw_sample_sheet_file_path) {
        error "[ERROR]: ❌ The --raw_sample_sheet_file_path parameter is required."
    } else {
        // Validate that the raw sample sheet file exists, is a file (not a directory), and has a .csv extension
        def raw_ss_path = file(params.raw_sample_sheet_file_path)
        if (!file(params.raw_sample_sheet_file_path).exists()) {
            error "[ERROR]: ❌ Raw sample sheet file does not exist: ${params.raw_sample_sheet_file_path}"
        }
        if (!raw_ss_path.isFile()) {
            error "[ERROR]: ❌ --raw_sample_sheet_file_path must be a file, not a directory: ${params.raw_sample_sheet_file_path}"
        }
        if (!raw_ss_path.name.toLowerCase().endsWith('.csv')) {
            error "[ERROR]: ❌ Raw sample sheet file must be a CSV file: ${params.raw_sample_sheet_file_path}"
        }
    }

    // Check for BCL directory
    if (!params.bcl_dir) {
        log.warn "[WARNING]: ⚠️ The --bcl_dir parameter is not specified. Default bcl_dir path will be set to: ${params.bcl_dir}"
    } else {
        // Validate that the BCL directory exists and is a directory
        def bcl_path = file(params.bcl_dir)
        if (!bcl_path.exists()) {
            error "[ERROR]: ❌ BCL directory does not exist: ${params.bcl_dir}"
        }
        if (!bcl_path.isDirectory()) {
            error "[ERROR]: ❌ --bcl_dir must be a directory, not a file: ${params.bcl_dir}"
        }
    }

    // Check for protocol prefix (required for batch naming)
    if (!params.protocol_prefix) {
        error "[ERROR]: ❌ The --protocol_prefix parameter is required for batch naming. Please provide a prefix that will be used to construct batch names, e.g., 'MIDAS2' or 'TecNante'."
    } else {
        // Validate that protocol_prefix is a non-empty string without spaces (to ensure valid batch names)
        if (params.protocol_prefix.trim().isEmpty()) {
            error "[ERROR]: ❌ The --protocol_prefix parameter cannot be empty. Please provide a valid prefix for batch naming, e.g., 'MIDAS2' or 'TecNante'."
        }
        if (params.protocol_prefix.contains(' ')) {
            error "[ERROR]: ❌ The --protocol_prefix parameter cannot contain spaces. Please provide a valid prefix for batch naming, e.g., 'MIDAS2' or 'TecNante'."
        }
    }

    // Check for batch IDs (required for identifying batches in the sample sheet and naming outputs)
    if (!params.batch_ids) {
        error "[ERROR]: ❌ The --batch_ids parameter is required for identifying batches in the sample sheet and naming outputs. Please provide a comma-separated list of batch IDs, e.g., '1,2,3'."
    } else {
        // Validate that batch_ids is a comma-separated list of integers
        def batch_ids_list = params.batch_ids.toString().split(',')
        if (batch_ids_list.size() == 0) {
            error "[ERROR]: ❌ The --batch_ids parameter must contain at least one batch ID. Please provide a comma-separated list of batch IDs, e.g., '1,2,3'."
        }
    }

    // Check genome (GEX) reference
    if (!params.path_ref_gex) {
        log.warn "[WARNING]: ⚠️ The --path_ref_gex parameter is not specified. Default GEX reference path will be set to: ${params.path_ref_gex}"
    } else {
        // Validate that the GEX reference path exists and is a directory
        def gex_ref_path = file(params.path_ref_gex)
        if (!gex_ref_path.exists()) {
            error "[ERROR]: ❌ GEX reference path does not exist: ${params.path_ref_gex}"
        }
        if (!gex_ref_path.isDirectory()) {
            error "[ERROR]: ❌ --path_ref_gex must be a directory, not a file: ${params.path_ref_gex}"
        }
    }

    // Check VDJ reference
    if (!params.path_ref_vdj) {
        error "[ERROR]: ❌ The --path_ref_vdj parameter is required for the VDJ reference used in Cellranger Multi. Please provide the path to the VDJ reference, e.g., '/path/to/refdata-vdj-GRCh38-alts-ensembl-2020-A'."
    } else {
        // Validate that the VDJ reference path exists and is a directory
        def vdj_ref_path = file(params.path_ref_vdj)
        if (!vdj_ref_path.exists()) {
            error "[ERROR]: ❌ VDJ reference path does not exist: ${params.path_ref_vdj}"
        }
        if (!vdj_ref_path.isDirectory()) {
            error "[ERROR]: ❌ --path_ref_vdj must be a directory, not a file: ${params.path_ref_vdj}"
        }
    }

    // Check output directory (optional, but if provided, must be a directory)
    if (!params.output_dir) {
        log.warn "[WARNING]: ⚠️ The --output_dir parameter is not specified. Default output directory will be set to: ${params.output_dir}"
    } else {
        def output_path = file(params.output_dir)
        if (!output_path.exists()) {
            log.warn "[WARNING]: ⚠️ Output directory does not exist. It will be created: ${params.output_dir}"
            output_path.mkdirs()
        } else {
            if (!output_path.isDirectory()) {
                error "[ERROR]: ❌ --output_dir must be a directory, not a file: ${params.output_dir}"
            }
        }
    }

    // Check preprocessing output directory (optional, but if provided, must be a directory)
    if (!params.preprocessing_output_dir) {
        log.warn "[WARNING]: ⚠️ The --preprocessing_output_dir parameter is not specified. Default preprocessing output directory will be set to: ${params.preprocessing_output_dir}"
    } else {
        def preproc_output_path = file(params.preprocessing_output_dir)
        if (!preproc_output_path.exists()) {
            log.warn "[WARNING]: ⚠️ Preprocessing output directory does not exist. It will be created: ${params.preprocessing_output_dir}"
            preproc_output_path.mkdirs()
        } else {
            if (!preproc_output_path.isDirectory()) {
                error "[ERROR]: ❌ --preprocessing_output_dir must be a directory, not a file: ${params.preprocessing_output_dir}"
            }
        }
    }

    // Check qc output directory (optional, but if provided, must be a directory)
    if (!params.qc_output_dir) {
        log.warn "[WARNING]: ⚠️ The --qc_output_dir parameter is not specified. Default QC output directory will be set to: ${params.qc_output_dir}"
    } else {
        def qc_output_path = file(params.qc_output_dir)
        if (!qc_output_path.exists()) {
            log.warn "[WARNING]: ⚠️ QC output directory does not exist. It will be created: ${params.qc_output_dir}"
            qc_output_path.mkdirs()
        } else {
            if (!qc_output_path.isDirectory()) {
                error "[ERROR]: ❌ --qc_output_dir must be a directory, not a file: ${params.qc_output_dir}"
            }
        }
    }

    // Check alignment output directory (optional, but if provided, must be a directory)
    if (!params.alignment_output_dir) {
        log.warn "[WARNING]: ⚠️ The --alignment_output_dir parameter is not specified. Default alignment output directory will be set to: ${params.alignment_output_dir}"
    } else {
        def align_output_path = file(params.alignment_output_dir)
        if (!align_output_path.exists()) {
            log.warn "[WARNING]: ⚠️ Alignment output directory does not exist. It will be created: ${params.alignment_output_dir}"
            align_output_path.mkdirs()
        } else {
            if (!align_output_path.isDirectory()) {
                error "[ERROR]: ❌ --alignment_output_dir must be a directory, not a file: ${params.alignment_output_dir}"
            }
        }
    }

    // Check MultiQC output directory (optional, but if provided, must be a directory)
    if (!params.multiqc_output_dir) {
        log.warn "[WARNING]: ⚠️ The --multiqc_output_dir parameter is not specified. Default MultiQC output directory will be set to: ${params.multiqc_output_dir}"
    } else {
        def multiqc_output_path = file(params.multiqc_output_dir)
        if (!multiqc_output_path.exists()) {
            log.warn "[WARNING]: ⚠️ MultiQC output directory does not exist. It will be created: ${params.multiqc_output_dir}"
            multiqc_output_path.mkdirs()
        } else {
            if (!multiqc_output_path.isDirectory()) {
                error "[ERROR]: ❌ --multiqc_output_dir must be a directory, not a file: ${params.multiqc_output_dir}"
            }
        }
    }

    // Check log directory (optional, but if provided, must be a directory)
    if (!params.log_dir) {
        log.warn "[WARNING]: ⚠️ The --log_dir parameter is not specified. Default log directory will be set to: ${params.log_dir}"
    } else {
        def log_path = file(params.log_dir)
        if (!log_path.exists()) {
            log.warn "[WARNING]: ⚠️ Log directory does not exist. It will be created: ${params.log_dir}"
            log_path.mkdirs()
        } else {
            if (!log_path.isDirectory()) {
                error "[ERROR]: ❌ --log_dir must be a directory, not a file: ${params.log_dir}"
            }
        }
    }

    // Check preprocessing log directory (optional, but if provided, must be a directory)
    if (!params.preprocessing_log_dir) {
        log.warn "[WARNING]: ⚠️ The --preprocessing_log_dir parameter is not specified. Default preprocessing log directory will be set to: ${params.preprocessing_log_dir}"
    } else {
        def preproc_log_path = file(params.preprocessing_log_dir)
        if (!preproc_log_path.exists()) {
            log.warn "[WARNING]: ⚠️ Preprocessing log directory does not exist. It will be created: ${params.preprocessing_log_dir}"
            preproc_log_path.mkdirs()
        } else {
            if (!preproc_log_path.isDirectory()) {
                error "[ERROR]: ❌ --preprocessing_log_dir must be a directory, not a file: ${params.preprocessing_log_dir}"
            }
        }
    }


    // Check QC log directory (optional, but if provided, must be a directory)
    if (!params.qc_log_dir) {        
        log.warn "[WARNING]: ⚠️ The --qc_log_dir parameter is not specified. Default QC log directory will be set to: ${params.qc_log_dir}"
    } else {
        def qc_log_path = file(params.qc_log_dir)
        if (!qc_log_path.exists()) {
            log.warn "[WARNING]: ⚠️ QC log directory does not exist. It will be created: ${params.qc_log_dir}"
            qc_log_path.mkdirs()
        } else {
            if (!qc_log_path.isDirectory()) {
                error "[ERROR]: ❌ --qc_log_dir must be a directory, not a file: ${params.qc_log_dir}"
            }
        }
    }

    // Check alignment log directory (optional, but if provided, must be a directory)    
    if (!params.alignment_log_dir) {
        log.warn "[WARNING]: ⚠️ The --alignment_log_dir parameter is not specified. Default alignment log directory will be set to: ${params.alignment_log_dir}"
    } else {
        def alignment_log_path = file(params.alignment_log_dir)
        if (!alignment_log_path.exists()) {
            log.warn "[WARNING]: ⚠️ Alignment log directory does not exist. It will be created: ${params.alignment_log_dir}"
            alignment_log_path.mkdirs()
        } else {
            if (!alignment_log_path.isDirectory()) {
                error "[ERROR]: ❌ --alignment_log_dir must be a directory, not a file: ${params.alignment_log_dir}"
            }
        }
    }

    // Check MultiQC log directory (optional, but if provided, must be a directory)
    if (!params.multiqc_log_dir) {
        log.warn "[WARNING]: ⚠️ The --multiqc_log_dir parameter is not specified. Default MultiQC log directory will be set to: ${params.multiqc_log_dir}"
    } else {
        def multiqc_log_path = file(params.multiqc_log_dir)
        if (!multiqc_log_path.exists()) {
            log.warn "[WARNING]: ⚠️ MultiQC log directory does not exist. It will be created: ${params.multiqc_log_dir}"
            multiqc_log_path.mkdirs()
        } else {
            if (!multiqc_log_path.isDirectory()) {
                error "[ERROR]: ❌ --multiqc_log_dir must be a directory, not a file: ${params.multiqc_log_dir}"
            }
        }
     }

    // Check localcores and/or localmemory (optional, but if provided, must be positive integers)
    if (!params.cpu_limit) {
        log.warn "[WARNING]: ⚠️ The --cpu_limit parameter is not specified. Default CPU limit will be set to: ${params.cpu_limit}"
    } else {
        if (params.cpu_limit <= 0) {
            error "[ERROR]: ❌ The --cpu_limit parameter must be a positive integer. Invalid value: ${params.cpu_limit}"
        } else {
            // If cpu_limit is specified, it should not exceed default cpu_limit (e.g., 16) to prevent overloading the system. This is a safeguard, but can be adjusted based on the specific environment and needs.
            if (params.cpu_limit > 16) {
                log.warn "[WARNING]: ⚠️ The specified --cpu_limit (${params.cpu_limit}) exceeds the recommended maximum of 16."
            }

        }
    }
    if (!params.memory_limit) {
        log.warn "[WARNING]: ⚠️ The --memory_limit parameter is not specified. Default memory limit will be set to: ${params.memory_limit}"
    } else {
        if (params.memory_limit <= 0) {
            error "[ERROR]: ❌ The --memory_limit parameter must be a positive integer. Invalid value: ${params.memory_limit}"
        } else {
            // If memory_limit is specified, it should not exceed default memory_limit (e.g., 64 GB) to prevent overloading the system. This is a safeguard, but can be adjusted based on the specific environment and needs.
            if (params.memory_limit > 64) {
                log.warn "[WARNING]: ⚠️ The specified --memory_limit (${params.memory_limit} GB) exceeds the recommended maximum of 64 GB."
            }
        }
    }

    // Log the successful completion of parameter validation
    log.info "[INFO]: ✅ Parameter validation passed."
}

// ========================================================================================
// MAIN WORKFLOW
// ========================================================================================

workflow {
    // -----------------------------------------------------------------------
    // Initialization and logging
    // -----------------------------------------------------------------------
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

    // Run PREPROCESSING module
    PREPROCESSING(
        file(params.raw_sample_sheet_file_path, checkIfExists: true), // Access the raw sample sheet path from params
        params.run_id, // Pass run_id for logging and output naming
        params.preprocessing_output_dir, // Pass output directory for standardized sample sheet
        params.preprocessing_log_dir, // Pass log directory for preprocessing logs
        params.today_date // Pass today's date for logging and output naming
    )

    // Capture outputs from PREPROCESSING
    ch_preprocessed_sample_sheet = PREPROCESSING.out.preprocessed_sample_sheet // Capture channel for standardized sample sheet
    ch_versions = ch_versions.mix(PREPROCESSING.out.versions) // Capture channel for versions information (mix for cumulation across steps)

    // -----------------------------------------------------------------------
    // STEP 2: Process BCL files to FASTQ according standardized sample sheet
    // -----------------------------------------------------------------------

    // Build mkfastq tuple input: (run_id, bcl_dir, preprocessed_sample_sheet)
    ch_mkfastq_input = ch_preprocessed_sample_sheet.map { preprocessed_sample_sheet ->
        tuple(params.run_id, file(params.bcl_dir), preprocessed_sample_sheet)
    }

    // Run CELLRANGER_MKFASTQ module
    CELLRANGER_MKFASTQ(
        ch_mkfastq_input, // Tuple input expected by module
        params.qc_output_dir, // Pass output directory for FASTQ files
        params.qc_log_dir, // Pass log directory for QC logs
        params.today_date // Pass today's date for logging and output naming
    )

    // Capture outputs from CELLRANGER_MKFASTQ
    ch_fastqs = CELLRANGER_MKFASTQ.out.fastqs // Capture channel for generated FASTQ files
    ch_versions = ch_versions.mix(CELLRANGER_MKFASTQ.out.versions) // Capture channel for versions information (mix for cumulation across steps)
    
    // -----------------------------------------------------------------------
    // STEP 3: Perform Alignment with Cellranger Multi
    // -----------------------------------------------------------------------

    // Convert the comma-separated batch IDs string into a list of batch id
    batch_ids_list = params.batch_ids.toString().split(',') // Split the comma-separated batch IDs into a list of batch_id

    // Create a channel with several elements (batch_id), each will be pass to a separate instance of CELLRANGER_MULTI, enabling parallel alignment across all requested batches.
    ch_batch_id = channel
        .from(batch_ids_list) // From the list of batch IDs
        .map { batch_id -> "${params.protocol_prefix}_batch${batch_id}" } // Map each batch_id in batch_ids_list to a full_batch_name

    // Build multi tuple input: (run_id, batch_id)
    ch_multi_input = ch_batch_id.map { batch_id -> tuple(params.run_id, batch_id) }
    
    // Run CELLRANGER_MULTI module
    CELLRANGER_MULTI(
        ch_multi_input, // Tuple input expected by module
        ch_fastqs, // Pass the channel of FASTQ directories generated from the previous step
        file(params.path_ref_gex, checkIfExists: true), // GEX reference from params
        file(params.path_ref_vdj, checkIfExists: true), // VDJ reference from params
        params.alignment_output_dir, // Pass output directory for Cellranger Multi results
        params.alignment_log_dir, // Pass log directory for Cellranger Multi
        params.today_date // Pass today's date for logging and output naming
    )

    // Capture outputs from CELLRANGER_MULTI
    ch_metrics = CELLRANGER_MULTI.out.metrics // Capture channel for metrics summary
    ch_web_summaries = CELLRANGER_MULTI.out.web_summaries // Capture channel for web summaries
    ch_versions = ch_versions.mix(CELLRANGER_MULTI.out.versions) // Capture channel for versions information (mix for cumulation across steps)

    // -----------------------------------------------------------------------
    // STEP 4: MultiQC report generation
    // -----------------------------------------------------------------------

    // Combine metrics summaries and web summaries into a single channel for MultiQC input
    ch_qc_files = 
        ch_metrics.map { _batch_id, f -> f } // In ch_metrics (that is a tuple(batch_id, metrics_file_path)), extract the file path (f) for each batch_id and create a channel of metrics summary files
        .mix(ch_web_summaries.map { _batch_id, f -> f }) 
        // In ch_web_summaries (that is a tuple(batch_id, web_summary_file_path)), 
        // extract the file path (f) for each batch_id and create a channel of web summary files 
        // ==> then mix it with the channel of metrics summary files to create a single channel of QC files for MultiQC input, containing both metrics summaries and web summaries from all batches

    // Run MULTIQC module
    MULTIQC(
        ch_qc_files.collect(), // Pass the channel of QC files collected (metrics summaries and web summaries) generated from Cellranger Multi
        params.run_id, // Pass run_id for logging and output naming
        params.multiqc_output_dir, // Pass output directory for MultiQC results
        params.multiqc_log_dir, // Pass log directory for MultiQC logs
        params.today_date // Pass today's date for logging and output naming
    )

    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    // -----------------------------------------------------------------------
    // STEP 5: Merge versions
    // -----------------------------------------------------------------------

    // Merge all per-module versions.yml files into one dated versions file
    MERGE_VERSIONS(
        params.run_id, // Pass run_id for logging and output naming
        ch_versions.collect(),
        params.run_traceability_log_dir, // Pass log directory for traceability logs (merged versions.yml)
        params.today_date
    )

    // -----------------------------------------------------------------------
    // Finishing workflow and logging
    // -----------------------------------------------------------------------

    // Capture values now because params can be unavailable in event handler scope
    def finalLogDir = params.log_dir
    def finalTodayDate = params.today_date
    def nfLogPath = '.nextflow.log'

    workflow.onComplete {
    def src = file(nfLogPath)
    def dst = file("${finalLogDir}/${finalTodayDate}_nextflow.log")

    if (src.exists()) {
        file(finalLogDir).mkdirs()
        src.copyTo(dst, overwrite: true)
        log.info "[INFO]: 📁 Nextflow log copied to: ${dst}"
    } else {
        log.warn "[WARNING]: ⚠️ Nextflow log not found: ${src}"
    }

    log.info "[INFO]: ✅✅✅ Pipeline completed successfully !"
    }

    workflow.onError {
        log.error "[ERROR]: ❌ Pipeline failed. Check .nextflow.log for details."
    } 

}
