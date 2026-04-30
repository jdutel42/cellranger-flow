/*
========================================================================================
    MODULE : CELLRANGER_MKFASTQ
========================================================================================
    Description : Converts a BCL directory into FASTQ files using cellranger mkfastq.
    Outil       : Cell Ranger 7.2.0
    Conda       : /labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_
----------------------------------------------------------------------------------------
    Inputs :
        tuple val(run_id), path(bcl_dir), path(sample_sheet)
    Outputs :
        path("fastq_output_${run_id}")          → ch_fastqs
        path "2_cellranger_mkfastq_versions.yml" → ch_versions
========================================================================================
*/

process CELLRANGER_MKFASTQ {

    // Tag and label for logging and resource allocation
    tag "CrF_${run_id}" // Log the run ID for traceability for "CrF" = "Cell Ranger mkFastq"
    label 'process_high' // Use a high resource label since this step can be computationally intensive

    // Container specification for reproducibility
    //container 'nfcore/cellranger:7.2.0'
    // Use conda environment for reproducibility
    conda '/labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_'

    // Publish FASTQ files to the output directory for downstream analysis
    publishDir (
        path    : params.qc_output_dir, // Use the QC output directory for FASTQ outputs
        mode    : 'copy', // Copy files to the output directory
        pattern : "fastq_output_${run_id}/**/*.fastq.gz", // Publish all files generated in work/fastq_output/ to qc_output_dir directory. ** is used to include all subdirectories (e.g., lane1, lane2) where FASTQ files are generated
    )
    
    // Publish logs to a dedicated directory
    publishDir (
        path    : params.qc_log_dir, // Use the QC log directory for logs
        mode    : 'copy', // Copy logs to the output directory
        pattern : "logs/*.log" // Publish log files generated in work/logs/ to qc_log_dir directory
    )

    // Declare process inputs
    input:
        tuple val(run_id), path(bcl_dir), path(preprocessed_sample_sheet)
        val qc_output_dir // Output directory for FASTQ files (val because it's STRING path)
        val qc_log_dir // Log directory for QC logs
        val today_date // Today's date for logging and output naming

    // Output the generated FASTQ files, versions information, and logs
    output:
        path("fastq_output_${run_id}"), emit: fastqs
        path "2_cellranger_mkfastq_versions.yml", emit: versions

    // Script section to run cellranger mkfastq
    script:
        """
        # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)
        set -euo pipefail

        # Create logs directory if it doesn't exist
        mkdir -p fastq_output_${run_id} logs

        # Redirect all log (stdout and stderr) to a log file for this process
        exec > >(tee -a logs/${today_date}_mkfastq_${run_id}.log) 2>&1

        # Load shared logging helpers
        source "${params.logging_script}"

        log_init "Step 2: Processing BCL files to FASTQ with Cellranger mkfastq for run_id = ${run_id}"
        log_log "Logs will be saved to ${qc_log_dir}/${today_date}_mkfastq_${run_id}.log"

        log_info "
        ╔═══════════════════════════════════════════════════════════════════════════════╗
        ║                         Cell Ranger mkfastq Process Script                    ║
        ╠═══════════════════════════════════════════════════════════════════════════════╣
        ║ Logging input parameters:
        ║ - run_id: ${run_id}
        ║ - bcl_dir: ${bcl_dir}
        ║ - sample_sheet: ${preprocessed_sample_sheet}
        ║ - cpus: ${params.cpu_limit}
        ║ - mem_gb: ${params.memory_limit}
        ║ - qc_output_dir: ${qc_output_dir}
        ║ - qc_log_dir: ${qc_log_dir}
        ║ - today_date: ${today_date}
        ╚═══════════════════════════════════════════════════════════════════════════════╝
        "

        # ============================================================================
        # Initialisation & Verification
        # ============================================================================

        log_verify "Verifying input files and directories..."

        # Verify that the BCL directory, preprocessed sample sheet directory exist before running cellranger mkfastq
        # If not, log an error message and exit with a non-zero status code to indicate failure
        if [ ! -d "${bcl_dir}" ]; then
            log_error "Missing BCL dir: ${bcl_dir}"
        fi

        if [ ! -f "${preprocessed_sample_sheet}" ]; then
            log_error "Missing sample sheet: ${preprocessed_sample_sheet}" >&2
        fi

        # Verify that the output directory exist, if not create it
        if [ ! -d "${qc_output_dir}" ]; then
            log_warning "Output directory ${qc_output_dir} does not exist. Creating it."
            mkdir -p "${qc_output_dir}"
        fi

        log_ok "Input verification completed successfully"

        # =============================================================================
        # Run Cellranger mkfastq
        # =============================================================================

        log_start "Running cellranger mkfastq for run_id = ${run_id}..."

        # See if run_ID is FlowCellID HCHNTDMX2 or CelineID 260323_A01789_0447_AHCHNTDMX2 (I think FlowCellID HCHNTDMX2 is better)
        /labos/UGM/dev/${params.cellranger_version}/bin/cellranger mkfastq \\
            --run="${bcl_dir}" \\
            --id="${run_id}" \\
            --csv="${preprocessed_sample_sheet}" \\
            --output-dir=fastq_output_${run_id} \\
            --localcores=${params.cpu_limit} \\
            --localmem=${params.memory_limit}

        log_ok "Cellranger mkfastq finished for run_id = ${run_id}"

        # -----------------------------------------------------------------------
        # Verify successful execution
        # -----------------------------------------------------------------------
        log_verify "Verifying cellranger mkfastq output..."

        if [ ! -d "fastq_output_${run_id}" ]; then
            log_error "Cellranger mkfastq did not generate the expected output directory: fastq_output_${run_id}/" >&2
        fi

        if ! find "fastq_output_${run_id}" -name "*.fastq.gz" | grep -q .; then
            log_error "No FASTQ files were generated in fastq_output_${run_id}/" >&2
        fi

        # Ensure the expected output directory exists and contains FASTQ files
        # Ensure at least one FASTQ file was generated
        FASTQ_COUNT=\$(find fastq_output_${run_id} -name "*.fastq.gz" | wc -l)
        if [ "\$FASTQ_COUNT" -eq 0 ]; then
            log_error "No FASTQ files generated in fastq_output_${run_id}/" >&2
        fi

        log_ok "Verified cellranger mkfastq output successfully"

        log_save "Cellranger mkfastq completed successfully for run_id = ${run_id} with \$FASTQ_COUNT FASTQ files generated in ${qc_output_dir}/fastq_output_${run_id}/"

        # -----------------------------------------------------------------------
        # Record tool version
        # -----------------------------------------------------------------------
        log_start "Recording tool versions for reproducibility..."

        cat <<-END_VERSIONS > 2_cellranger_mkfastq_versions.yml
        "${task.process}":
            cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
        END_VERSIONS

        log_ok "Tool versions recorded successfully in 2_cellranger_mkfastq_versions.yml"

        # -----------------------------------------------------------------------
        # End
        # -----------------------------------------------------------------------
        
        log_save "Cellranger multi output fo run_id ${run_id} saved to ${qc_output_dir}/fastq_output_${run_id}/."
        log_log "Versions information will be saved to ${params.run_traceability_log_dir}/${today_date}_versions.yaml"
        log_log "Logs saved to ${qc_log_dir}/${today_date}_mkfastq_${run_id}.log"

        log_success "Cell Ranger mkfastq processed BCL to FASTQ successfully for run_id = ${run_id} and FASTQ files are available at ${params.qc_output_dir} !"
        """
    
    stub:
    """
    mkdir -p fastq_output_${run_id}/lane1 logs
    
    # Create dummy FASTQ files
    touch fastq_output_${run_id}/lane1/${run_id}_S1_L001_R1_001.fastq.gz
    touch fastq_output_${run_id}/lane1/${run_id}_S1_L001_R2_001.fastq.gz
    
    # Create minimal versions file
    cat > 2_cellranger_mkfastq_versions.yml <<EOF
"CELLRANGER_MKFASTQ":
    "cellranger": "${params.cellranger_version}"
    "run_id": "${run_id}"
EOF
    """
}
