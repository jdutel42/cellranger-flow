/*
========================================================================================
    SUB-WORKFLOW : GEX_VDJ_subworkflow
========================================================================================
    Description : Orchestrates the sample sheet preprocessing, generation and processingof FASTQ files for GEX and VDJ data using CellRanger Mkfastq and Multi, 
                  followed by QC report generation with MultiQC and merging of versions information.
    Steps       :
        1. PREPROCESSING       (raw sample sheet -> standardized sample sheet)
        2. CELLRANGER_MKFASTQ  (BCL + standardized sample sheet -> FASTQ per sample)
        3. CELLRANGER_MULTI    (FASTQ per sample + references -> matrices + metrics + web summaries)
        4. MULTIQC             (metrics summaries + web summaries -> MultiQC report)
        5. MERGE_VERSIONS      (versions.yml from all steps -> merged versions.yml with run_id and date)
========================================================================================
*/

// Activate DSL2 language, which is required for modular workflows and subworkflows
nextflow.enable.dsl = 2

// ========================================================================================
// IMPORTS MODULS
// ========================================================================================

include { PREPROCESSING       } from '../modules/preprocessing.nf'
include { CELLRANGER_MKFASTQ  } from '../modules/cellranger_mkfastq.nf'
include { CELLRANGER_MULTI    } from '../modules/cellranger_multi.nf'
include { MULTIQC             } from '../modules/multiqc.nf'
include { MERGE_VERSIONS      } from '../modules/merge_versions.nf'

workflow FASTQ_TO_COUNT {

    take:
        ch_fastq_dirs   // Channel[ tuple(sample_id, fastq_dir) ]

    main:
    // -----------------------------------------------------------------------
    // Initialization and logging
    // -----------------------------------------------------------------------
    // Print the date and time when the workflow starts
    log.info "[INFO]: 🚀🚀🚀 Starting CellRanger_Flow workflow at ${params.today_date} with run_id ${params.run_id}..."

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
    if (shouldRun("preprocessing")) {
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
    } else {
        ch_preprocessed_sample_sheet = channel.fromPath(params.preprocessed_sample_sheet_path)
    }

    // -----------------------------------------------------------------------
    // STEP 2: Process BCL files to FASTQ according standardized sample sheet
    // -----------------------------------------------------------------------

    if (shouldRun("mkfastq")) {
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
    } else {
        ch_fastqs = channel.fromPath(params.fastqs_dir_path)
    }

    // -----------------------------------------------------------------------
    // STEP 3: Perform Alignment with Cellranger Multi
    // -----------------------------------------------------------------------

    if (shouldRun("multi")) {
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
            params.genome_reference_path, // GEX reference from params
            params.vdj_reference_path, // VDJ reference from params
            params.alignment_output_dir, // Pass output directory for Cellranger Multi results
            params.alignment_log_dir, // Pass log directory for Cellranger Multi
            params.today_date // Pass today's date for logging and output naming
        )

        // Capture outputs from CELLRANGER_MULTI
        ch_metrics = CELLRANGER_MULTI.out.metrics // Capture channel for metrics summary
        ch_web_summaries = CELLRANGER_MULTI.out.web_summaries // Capture channel for web summaries
        ch_versions = ch_versions.mix(CELLRANGER_MULTI.out.versions) // Capture channel for versions information (mix for cumulation across steps)
    } else {
        ch_metrics = channel.fromPath(params.metrics_dir_path)
        ch_web_summaries = channel.fromPath(params.web_summaries_dir_path)
    }
    // -----------------------------------------------------------------------
    // STEP 4: MultiQC report generation
    // -----------------------------------------------------------------------

    if (shouldRun("multiqc")) {
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
            ch_fastqs.collect(), // Pass the channel of FASTQ directories generated from the previous step for MultiQC to link raw data in the report
            params.run_id, // Pass run_id for logging and output naming
            params.multiqc_output_dir, // Pass output directory for MultiQC results
            params.multiqc_log_dir, // Pass log directory for MultiQC logs
            params.today_date // Pass today's date for logging and output naming
        )

        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    } else {
        // If MultiQC step is skipped, we can still capture the versions information from previous steps
        ch_versions = ch_versions.mix(channel.fromPath(params.versions_dir_path))
    }

    // -----------------------------------------------------------------------
    // STEP 5: Merge versions
    // -----------------------------------------------------------------------

    if (shouldRun("merge_versions")) {
        // Merge all per-module versions.yml files into one dated versions file
        MERGE_VERSIONS(
            params.run_id, // Pass run_id for logging and output naming
            ch_versions.collect(),
            params.run_traceability_log_dir, // Pass log directory for traceability logs (merged versions.yml)
            params.today_date
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

    emit:
        matrices      = CELLRANGER_COUNT.out.matrices
        metrics       = CELLRANGER_COUNT.out.metrics
        web_summaries = CELLRANGER_COUNT.out.web_summaries
        molecule_info = CELLRANGER_COUNT.out.molecule_info
        versions      = ch_versions
}
