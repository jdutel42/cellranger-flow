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

include { PROCESSING_SAMPLE_SHEET } from './modules/processing_sample_sheet.nf'
include { CELLRANGER_MKFASTQ  } from './modules/cellranger_mkfastq.nf'
include { CELLRANGER_MULTI    } from './modules/cellranger_multi.nf'
include { MULTIQC             } from './modules/multiqc.nf'
include { MERGE_VERSIONS      } from './modules/merge_versions.nf'

// ========================================================================================
// IMPORTS SUB-WORKFLOWS
// ========================================================================================

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
    ║  protocol_id              : ${params.protocol_id}
    ║  batch_ids                : ${params.batch_ids}
    ║  species                  : ${params.species}
    ║  genome_version           : ${params.genome_version}
    ║  genome_reference_path    : ${params.genome_reference_path}
    ║  vdj_reference_path       : ${params.vdj_reference_path}
    ║  output_dir               : ${params.output_dir}
    ║  processing_output_dir    : ${params.processing_output_dir}
    ║  qc_output_dir            : ${params.qc_output_dir}
    ║  alignment_output_dir     : ${params.alignment_output_dir}
    ║  multiqc_output_dir       : ${params.multiqc_output_dir}
    ║  log_dir                  : ${params.log_dir}
    ║  cpu/memory_limit         : ${params.cpu_limit} / ${params.memory_limit}GB
    ║  pipeline_version         : ${params.pipeline_version}
    ║  min_nextflow_version     : ${params.min_nextflow_version}
    ║  cellranger_version       : ${params.cellranger_version}
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

    // Check for protocol ID (required for batch naming)
    if (!params.protocol_id) {
        error "[ERROR]: ❌ The --protocol_id parameter is required for batch naming. Please provide an ID that will be used to construct batch names, e.g., 'MIDAS2' or 'TecNante'."
    } else {
        // Validate that protocol_id is a non-empty string without spaces (to ensure valid batch names)
        if (params.protocol_id.trim().isEmpty()) {
            error "[ERROR]: ❌ The --protocol_id parameter cannot be empty. Please provide a valid ID for batch naming, e.g., 'MIDAS2' or 'TecNante'."
        }
        if (params.protocol_id.contains(' ')) {
            error "[ERROR]: ❌ The --protocol_id parameter cannot contain spaces. Please provide a valid ID for batch naming, e.g., 'MIDAS2' or 'TecNante'."
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

    // Check species
    if (!params.species) {
        error "[ERROR]: ❌ The --species parameter is not specified. Please provide a valid species, either 'human' or 'mouse'."
    } else {
        // Validate that the specified species is one of the allowed values
        if (!(params.species in ["human", "mouse"])) {
            error "[ERROR]: ❌ The --species parameter must be either 'human' or 'mouse'. Invalid value: ${params.species}"
        }
    }
    if (params.species == "human" && !params.genome_reference_path) {
        log.warn "[WARNING]: ⚠️ The --genome_reference_path parameter is not specified. Default and latest human GEX reference path will be set to: ${params.genome_reference_path}"
    }
    if (params.species == "mouse" && !params.genome_reference_path) {
        log.warn "[WARNING]: ⚠️ The --genome_reference_path parameter is not specified. Default mouse GEX reference path will be set to: ${params.genome_reference_path}"
    }
    

    // Check genome (GEX) reference
    if (!params.genome_reference_path) {
        log.warn "[WARNING]: ⚠️ The --genome_reference_path parameter is not specified. Default GEX reference path will be set according to the specified species and genome version : ${params.species} ${params.genome_version}"
        if (!params.species) {
            error "[ERROR]: ❌ The --species parameter is required to set the default GEX reference path. Please provide a valid species, either 'human' or 'mouse'."
        } else if (params.species == "human" || params.species == "Homo sapiens") {
            if (!params.genome_version) {
                error "[ERROR]: ❌ The --genome_version parameter is required for the human species. Please provide a valid genome version, e.g., 'hg38' or 'GRCh38'."
            } else if (params.genome_version == "hg38" || params.genome_version == "GRCh38") {
                params.genome_reference_path = "/labos/UGM/dev/cellranger-pipe/refdata-gex-GRCh38-2020-A"
            } else if (params.genome_version == "hg19" || params.genome_version == "GRCh37") {
                params.genome_reference_path = "/labos/UGM/dev/cellranger-pipe/refdata-gex-GRCh37-2020-A"
            } else {
                error "[ERROR]: ❌ Unsupported genome version for human species: ${params.genome_version}. Supported version is: GRCh38."
            }
        } else if (params.species == "mouse" || params.species == "mus_musculus") {
            if (!params.genome_version) {
                error "[ERROR]: ❌ The --genome_version parameter is required for the mouse species. Please provide a valid genome version, e.g., 'mm39' or 'GRCm39'."
            } else if (params.genome_version == "mm39" || params.genome_version == "GRCm39") {
                params.genome_reference_path = "/home/dutel/Data/Reference/Other/GEX/refdata-gex-GRCm39-2024-A"
            } else {
                error "[ERROR]: ❌ Unsupported genome version for mouse species: ${params.genome_version}. Supported version is: mm39."
            }
        } else {
            error "[ERROR]: ❌ Unsupported species: ${params.species}. Supported species are: human and mouse."
        }
    } else {
        // Validate that the GEX reference path exists and is a directory
        def gex_ref_path = file(params.genome_reference_path)
        if (!gex_ref_path.exists()) {
            error "[ERROR]: ❌ GEX reference path does not exist: ${params.genome_reference_path}"
        }
        if (!gex_ref_path.isDirectory()) {
            error "[ERROR]: ❌ --genome_reference_path must be a directory, not a file: ${params.genome_reference_path}"
        }
    }

    // Check VDJ reference
    if (!params.vdj_reference_path) {
        log.warn "[WARNING]: ⚠️ The --vdj_reference_path parameter is not specified. Default VDJ reference path will be set according to the specified species and genome version : ${params.species} ${params.genome_version}"
        if (!params.species) {
            error "[ERROR]: ❌ The --species parameter is required to set the default VDJ reference path. Please provide a valid species, either 'human' or 'mouse'."
        } else if (params.species == "human" || params.species == "Homo sapiens") {
            if (!params.genome_version) {
                error "[ERROR]: ❌ The --genome_version parameter is required for the human species. Please provide a valid genome version, e.g., 'hg38' or 'GRCh38'."
            } else if (params.genome_version == "hg38" || params.genome_version == "GRCh38") {
                params.vdj_reference_path = "/labos/UGM/dev/cellranger-pipe/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0"
            } else if (params.genome_version == "hg19" || params.genome_version == "GRCh37") {
                params.vdj_reference_path = "/labos/UGM/dev/cellranger-pipe/refdata-cellranger-vdj-GRCh37-alts-ensembl-2020-A"
            } else {
                error "[ERROR]: ❌ Unsupported genome version for human species: ${params.genome_version}. Supported version is: GRCh38."
            }
        } else if (params.species == "mouse" || params.species == "mus_musculus") {
            if (!params.genome_version) {
                error "[ERROR]: ❌ The --genome_version parameter is required for the mouse species. Please provide a valid genome version, e.g., 'mm39' or 'GRCm39'."
            } else if (params.genome_version == "mm39" || params.genome_version == "GRCm39") {
                params.vdj_reference_path = "/home/dutel/Data/Reference/Other/VDJ/refdata-vdj-GRCm39-2024-A"
            } else {
                error "[ERROR]: ❌ Unsupported genome version for mouse species: ${params.genome_version}. Supported version is: mm39."
            }
        } else {
            error "[ERROR]: ❌ Unsupported species: ${params.species}. Supported species are: human and mouse."
        }
    } else {
        // Validate that the VDJ reference path exists and is a directory
        def vdj_ref_path = file(params.vdj_reference_path)
        if (!vdj_ref_path.exists()) {
            error "[ERROR]: ❌ VDJ reference path does not exist: ${params.vdj_reference_path}"
        }
        if (!vdj_ref_path.isDirectory()) {
            error "[ERROR]: ❌ --vdj_reference_path must be a directory, not a file: ${params.vdj_reference_path}"
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

    // Check processing output directory (optional, but if provided, must be a directory)
    if (!params.processing_output_dir) {
        log.warn "[WARNING]: ⚠️ The --processing_output_dir parameter is not specified. Default processing output directory will be set to: ${params.processing_output_dir}"
    } else {
        def preproc_output_path = file(params.processing_output_dir)
        if (!preproc_output_path.exists()) {
            log.warn "[WARNING]: ⚠️ Processing output directory does not exist. It will be created: ${params.processing_output_dir}"
            preproc_output_path.mkdirs()
        } else {
            if (!preproc_output_path.isDirectory()) {
                error "[ERROR]: ❌ --processing_output_dir must be a directory, not a file: ${params.processing_output_dir}"
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
    if (!params.processing_log_dir) {
        log.warn "[WARNING]: ⚠️ The --processing_log_dir parameter is not specified. Default processing log directory will be set to: ${params.processing_log_dir}"
    } else {
        def preproc_log_path = file(params.processing_log_dir)
        if (!preproc_log_path.exists()) {
            log.warn "[WARNING]: ⚠️ Processing log directory does not exist. It will be created: ${params.processing_log_dir}"
            preproc_log_path.mkdirs()
        } else {
            if (!preproc_log_path.isDirectory()) {
                error "[ERROR]: ❌ --processing_log_dir must be a directory, not a file: ${params.processing_log_dir}"
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
            if (params.memory_limit > 80) {
                log.warn "[WARNING]: ⚠️ The specified --memory_limit (${params.memory_limit} GB) exceeds the recommended maximum of 64 GB."
            }
        }
    }

    if (!params.cellranger_version) {
        params.cellranger_version = 'cellranger-7.1.0'
        log.warn "[WARNING]: ⚠️ The --cellranger_version parameter is not specified. Default Cell Ranger version will be set to: ${params.cellranger_version}"
    } else {
        // Validate that the specified Cell Ranger version is one of the allowed values
        if (!(params.cellranger_version in ["cellranger-7.1.0", "cellranger-9.0.1"])) {
            error "[ERROR]: ❌ The --cellranger_version parameter must be cellranger-7.1.0 or cellranger-9.0.1. Invalid value: ${params.cellranger_version}"
        }
    }

    // Log the successful completion of parameter validation
    log.info "[INFO]: ✅ Parameter validation passed."
}



def shouldRun(String step) {
    // return params.mode == "full" || params.mode == step // Condensed version
    // OR
    if (params.mode == "full") {
        return true
    } else if (params.mode == step) {
        return true
    } else {
        return false
    }
}





// ========================================================================================
// MAIN WORKFLOW
// ========================================================================================

workflow {
    // -----------------------------------------------------------------------
    // Initialization and logging
    // -----------------------------------------------------------------------
    // Print the date and time when the workflow starts
    log.info "[INFO]: 🚀🚀🚀 Starting CellRanger_Flow workflow with run_id ${params.run_id}..."

    // Print header and validate parameters at the start of the workflow
    printHeader()

    // Validate parameters before starting any processing
    validateParams()

    // Channel for version files (for MultiQC and traceability)
    // Initialization is needed to ensure accumulation versions from multiple steps
    ch_versions = channel.empty()

    // -----------------------------------------------------------------------
    // STEP 1: Process the sample sheet to standardize it for downstream tools
    // -----------------------------------------------------------------------

    // Run PROCESSING_SAMPLE_SHEET module
    if (shouldRun("processing_sample_sheet")) {
        PROCESSING_SAMPLE_SHEET(
            file(params.raw_sample_sheet_file_path, checkIfExists: true) // Access the raw sample sheet path from params
        )
        
        // Capture outputs from PROCESSING_SAMPLE_SHEET
        ch_processed_sample_sheet = PROCESSING_SAMPLE_SHEET.out.ch_processed_sample_sheet // Capture channel for standardized sample sheet
        ch_versions = ch_versions.mix(PROCESSING_SAMPLE_SHEET.out.ch_versions) // Capture channel for versions information (mix for cumulation across steps)
    } else {
        ch_processed_sample_sheet = channel.fromPath(params.processed_sample_sheet_path)
    }

    // -----------------------------------------------------------------------
    // STEP 2: Process BCL files to FASTQ according standardized sample sheet
    // -----------------------------------------------------------------------

    // Run CELLRANGER_MKFASTQ module
    if (shouldRun("cellranger_mkfastq")) {
        CELLRANGER_MKFASTQ(
            ch_processed_sample_sheet // Pass the channel of processed sample sheet generated from the previous step
        )

        // Capture outputs from CELLRANGER_MKFASTQ
        ch_fastqs = CELLRANGER_MKFASTQ.out.ch_fastqs.map { file(it).toRealPath() } // Capture channel (absolute paths) for generated FASTQ files
        ch_versions = ch_versions.mix(CELLRANGER_MKFASTQ.out.ch_versions) // Capture channel for versions information (mix for cumulation across steps)
    } else {
        ch_fastqs = channel.fromPath(params.fastqs_dir_path)
    }

    // -----------------------------------------------------------------------
    // STEP 3: Perform Alignment with Cellranger Multi
    // -----------------------------------------------------------------------

    if (shouldRun("cellranger_multi")) {
        // Use sample_ids directly when provided; otherwise reconstruct batch names from batch_ids
        //def multi_sample_ids = params.sample_ids
        //    ? params.sample_ids.toString().split(',').collect { sample_id -> sample_id.trim() }
        //    : params.batch_ids.toString().split(',').collect { batch_id -> "${params.protocol_id}_batch${batch_id}" }

        // Build multi tuple input: (run_id, sample_id/batch_id)
        //ch_multi_input = channel.from(multi_sample_ids).map { sample_id -> tuple(params.run_id, sample_id) }
        
        // Create a channel of sample IDs for Cellranger Multi
        ch_sample_ids = Channel.from(params.sample_ids.split(',').collect { it.trim() })

        // Create channel of tuples (fastq_dir, sample_id) with combine operator
        ch_multi_input = ch_fastqs.combine(ch_sample_ids)

        // Run CELLRANGER_MULTI module
        CELLRANGER_MULTI(
            ch_multi_input // Tuples input expected by module
        )

        // Capture outputs from CELLRANGER_MULTI
        ch_multiqc_input = CELLRANGER_MULTI.out.ch_multi_output // Capture channel for MultiQC input (metrics summaries and web summaries)
        ch_metrics = CELLRANGER_MULTI.out.ch_metrics // Capture channel for metrics summary
        ch_web_summaries = CELLRANGER_MULTI.out.ch_web_summaries // Capture channel for web summaries
        ch_versions = ch_versions.mix(CELLRANGER_MULTI.out.ch_versions) // Capture channel for versions information (mix for cumulation across steps)
    } else {
        ch_metrics = channel.fromPath(params.metrics_dir_path)
        ch_web_summaries = channel.fromPath(params.web_summaries_dir_path)
    }

    // -----------------------------------------------------------------------
    // STEP 4: MultiQC report generation
    // -----------------------------------------------------------------------

    if (shouldRun("multiqc")) {
        // Combine metrics summaries and web summaries into a single channel for MultiQC input
        ch_multi_output = ch_metrics.map { _, file -> file }.mix(ch_web_summaries.map { _, file -> file }).mix(ch_fastqs)
        ch_multiqc_input = ch_multi_output.collect() // Collect all files into a single list for MultiQC input

        // Run MULTIQC module
        MULTIQC(
            ch_multiqc_input
        )

        ch_versions = ch_versions.mix(MULTIQC.out.ch_versions)
    } else {
        // If MultiQC step is skipped, we can still capture the versions information from previous steps
        ch_versions = ch_versions.mix(channel.fromPath(params.versions_dir_path))
    }

    // -----------------------------------------------------------------------
    // STEP 5: Merge versions
    // -----------------------------------------------------------------------

    if (shouldRun("merge_versions")) {
        // Merge all per-module versions.yml files into one dated versions file
        ch_versions = ch_versions.collect() // Collect all versions files into a single list for merging
        MERGE_VERSIONS(
            ch_versions
        )
    } else {
        // If merge versions step is skipped, we can still capture the versions information from previous steps
        ch_versions = ch_versions.mix(channel.fromPath(params.versions_dir_path))
    }

    // -----------------------------------------------------------------------
    // Finishing workflow and logging
    // -----------------------------------------------------------------------

    // Capture values now because params can be unavailable in event handler scope
    def finalLogDir = params.log_dir
    // def finalTodayDate = params.today_date
    def nfLogPath = '.nextflow.log'

    workflow.onComplete {
    def src = file(nfLogPath)
    def dst = file("${finalLogDir}/nextflow.log")

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