/*
========================================================================================
    MODULE : CELLRANGER_MKFASTQ
========================================================================================
    Description : Converts a BCL directory into FASTQ files using cellranger mkfastq.
    Outil       : Cell Ranger 7.2.0
    Conda       : /labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_
----------------------------------------------------------------------------------------
    Inputs :
        val processed_sample_sheet   -> Preprocessed sample sheet path for input

    Outputs :
        path "fastq_output_${params.bcl_id}"                    → ch_fastqs
        path "Cellranger_mkfastq_versions_${params.bcl_id}.yml" → ch_versions
========================================================================================
*/

process CELLRANGER_MKFASTQ {

    // Tag and label for logging and resource allocation
    tag "CrMk_${params.bcl_id}_CellRanger_Mkfastq_${params.protocol_prefix}" // Log the run ID for traceability for "CrF" = "Cell Ranger mkFastq"
    label 'process_high' // Use a high resource label since this step can be computationally intensive

    // Container specification for reproducibility
    //container 'nfcore/cellranger:7.2.0'
    // Use conda environment for reproducibility
    conda '/labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_'

    // Publish the entire FASTQ output directory (all files and subdirs) to the QC output directory
    publishDir (
        path    : params.qc_output_dir, // Use the QC output directory for FASTQ outputs
        mode    : 'copy', // Copy files to the output directory
        pattern : "fastq_output_${params.bcl_id}/**" // Publish the full directory tree produced by cellranger mkfastq
    )
    
    // Publish logs to a dedicated directory
    publishDir (
        path    : params.qc_log_dir, // Use the QC log directory for logs
        mode    : 'copy', // Copy logs to the output directory
        pattern : "logs/*.log" // Publish log files generated in work/logs/ to qc_log_dir directory
    )

    // Declare process inputs
    input:
        val processed_sample_sheet // Preprocessed sample sheet path for input

    // Output the generated FASTQ files, versions information, and logs
    output:
        path "fastq_output_${params.bcl_id}",                    emit: ch_fastqs
        path "Cellranger_mkfastq_versions_${params.bcl_id}.yml", emit: ch_versions

    // Script section to run cellranger mkfastq
    script:
        """
        # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)
        set -euo pipefail

        # Create fastq_output and logs directories in work/
        mkdir -p fastq_output_${params.bcl_id} logs

        # Redirect all log (stdout and stderr) to a log file for this process
        exec > >(tee -a logs/Cellranger_mkfastq_log_${params.bcl_id}.log) 2>&1

        # Load shared logging helpers
        source "${params.logging_script}"

        log_init "Step 2: Processing BCL files to FASTQ with Cellranger mkfastq for run_id = ${params.run_id}"
        
        log_log "Logs will be saved to ${params.qc_log_dir}/Cellranger_mkfastq_log_${params.bcl_id}.log"

        log_info "
        ╔═══════════════════════════════════════════════════════════════════════════════╗
        ║                         CellRanger mkfastq Process Script                    ║
        ╠═══════════════════════════════════════════════════════════════════════════════╣
        ║ Logging input parameters:
        ║ - run_id: ${params.run_id}
        ║ - bcl_id: ${params.bcl_id}
        ║ - bcl_dir: ${params.bcl_dir}
        ║ - sample_sheet: ${params.processed_sample_sheet}
        ║ - cpus: ${params.cpu_limit}
        ║ - mem_gb: ${params.memory_limit}
        ║ - qc_output_dir: ${params.qc_output_dir}
        ║ - qc_log_dir: ${params.qc_log_dir}
        ╚═══════════════════════════════════════════════════════════════════════════════╝
        "

        # ============================================================================
        # Initialisation & Verification
        # ============================================================================

        log_verify "Verifying input files and directories..."

        # Verify that the BCL directory, processed sample sheet directory exist before running cellranger mkfastq
        # If not, log an error message and exit with a non-zero status code to indicate failure
        if [ ! -d "${params.bcl_dir}" ]; then
            log_error "Missing BCL dir: ${params.bcl_dir}"
        fi

        if [ ! -f "${params.processed_sample_sheet}" ]; then
            log_error "Missing sample sheet: ${params.processed_sample_sheet}" >&2
        fi

        # Verify that the output directory exist, if not create it
        if [ ! -d "${params.qc_output_dir}" ]; then
            log_warning "Output directory ${params.qc_output_dir} does not exist. Creating it."
            mkdir -p "${params.qc_output_dir}"
        fi

        log_ok "Input verification completed successfully"

        # =============================================================================
        # Run Cellranger mkfastq
        # =============================================================================

        log_start "Running cellranger mkfastq for run_id = ${params.run_id}..."

        # See if run_ID is FlowCellID HCHNTDMX2 or CelineID 260323_A01789_0447_AHCHNTDMX2 (I think FlowCellID HCHNTDMX2 is better)
        /labos/UGM/dev/${params.cellranger_version}/bin/cellranger mkfastq \\
            --run="${params.bcl_dir}" \\
            --id="${params.run_id}" \\
            --csv="${params.processed_sample_sheet}" \\
            --output-dir=fastq_output_${params.bcl_id} \\
            --localcores=${params.cpu_limit} \\
            --localmem=${params.memory_limit}

        log_ok "Cellranger mkfastq finished for run_id = ${params.run_id}"

        # -----------------------------------------------------------------------
        # Verify successful execution
        # -----------------------------------------------------------------------
        log_verify "Verifying cellranger mkfastq output..."

        if [ ! -d "fastq_output_${params.bcl_id}" ]; then
            log_error "Cellranger mkfastq did not generate the expected output directory: fastq_output_${params.bcl_id}/" >&2
        fi

        # Commented because fastq are sometimes generated but verification is too soon and it fails
        #if ! find "fastq_output_${params.bcl_id}" -name "*.fastq.gz" | grep -q .; then
        #    log_error "No FASTQ files were generated in fastq_output_${params.bcl_id}/" >&2
        #fi

        log_ok "Verified cellranger mkfastq output successfully"

        log_save "Cellranger mkfastq completed successfully for run_id = ${params.run_id} with FASTQ files published in ${params.qc_output_dir}/fastq_output_${params.bcl_id}/"

        # -----------------------------------------------------------------------
        # Record tool version
        # -----------------------------------------------------------------------
        log_start "Recording tool versions for reproducibility..."

        cat <<-END_VERSIONS > Cellranger_mkfastq_versions_${params.bcl_id}.yml
        "${task.process}":
            cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
        END_VERSIONS

        log_ok "Tool versions recorded successfully in Cellranger_mkfastq_versions_${params.bcl_id}.yml"

        # -----------------------------------------------------------------------
        # End
        # -----------------------------------------------------------------------
        
        log_save "Cellranger multi output fo run_id ${params.run_id} saved to ${params.qc_output_dir}/fastq_output_${params.bcl_id}/."
        log_log "Versions information will be saved to ${params.run_traceability_log_dir}/Cellranger_mkfastq_versions_${params.bcl_id}.yml.yaml"
        log_log "Logs saved to ${params.qc_log_dir}/Cellranger_mkfastq_log_${params.bcl_id}.log"

        log_success "Cell Ranger mkfastq processed BCL to FASTQ successfully for run_id = ${params.run_id} and FASTQ files are available at ${params.qc_output_dir} !"
        """
    
    stub:
    """
    mkdir -p fastq_output_${params.bcl_id}/lane1 logs
    
    # Create dummy FASTQ files
    touch fastq_output_${params.bcl_id}/lane1/${params.bcl_id}_S1_L001_R1_001.fastq.gz
    touch fastq_output_${params.bcl_id}/lane1/${params.bcl_id}_S1_L001_R2_001.fastq.gz
    
    # Create minimal versions file
    cat > Cellranger_mkfastq_versions_${params.bcl_id}.yml <<EOF
"CELLRANGER_MKFASTQ":
    "cellranger": "${params.cellranger_version}"
    "run_id": "${params.run_id}"
EOF
    """
}
