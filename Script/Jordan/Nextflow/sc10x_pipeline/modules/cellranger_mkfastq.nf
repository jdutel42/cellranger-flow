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
        tuple val(run_id), path("fastq_output") → ch_fastqs
        path "logs/*.log"                    → ch_logs
        path "versions.yml"                             → ch_versions
========================================================================================
*/

process CELLRANGER_MKFASTQ {

    // Tag and label for logging and resource allocation
    tag "mkfastq | run_id: ${run_id}" // Log the run ID for traceability
    label 'process_low' // Use a high resource label since this step can be computationally intensive

    // Container specification for reproducibility
    //container 'nfcore/cellranger:7.2.0'
    // Use conda environment for reproducibility
    conda '/labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_'

    // Publish FASTQ files to the output directory for downstream analysis
    publishDir (
        path    : qc_output_dir, // Use the QC output directory for FASTQ outputs
        mode    : 'copy', // Copy files to the output directory
        pattern : "fastq_output/**/*.fastq.gz", // Publish all files generated in work/fastq_output/ to qc_output_dir directory. ** is used to include all subdirectories (e.g., lane1, lane2) where FASTQ files are generated
    )
    
    // Publish logs to a dedicated directory
    publishDir (
        path    : qc_log_dir, // Use the QC log directory for logs
        mode    : 'copy', // Copy logs to the output directory
        pattern : "logs/*.log" // Publish all log files generated in work/logs/ to qc_log_dir directory
    )

    // Declare process inputs
    input:
        tuple val(run_id), path(bcl_dir), path(preprocessed_sample_sheet)
        val qc_output_dir // Output directory for FASTQ files (val because it's STRING path)
        val qc_log_dir // Log directory for QC logs

    // Output the generated FASTQ files, versions information, and logs
    output:
        tuple val(run_id), path("fastq_output"), emit: fastqs
        path "logs/*.log", emit: logs
        path "versions.yml", emit: versions

    // Script section to run cellranger mkfastq
    script:
    """
    # Stop the script if a command fails, if an undefined variable is used, or if a command in a pipeline fails
    set -euo pipefail 

    # Redirect all log (stdout and stderr) to a log file for this process
    mkdir -p fastq_output logs
    exec > >(tee -a logs/mkfastq_${run_id}.log) 2>&1

    ##########################################
    #                 Verif                  #
    ##########################################

    # Verify that the BCL directory, preprocessed sample sheet directory exist before running cellranger mkfastq
    # If not, log an error message and exit with a non-zero status code to indicate failure
    if [ ! -d "${bcl_dir}" ]; then
        echo "[ERROR] Missing BCL dir: ${bcl_dir}" >&2
        exit 1
    fi

    if [ ! -f "${preprocessed_sample_sheet}" ]; then
        echo "[ERROR] Missing sample sheet: ${preprocessed_sample_sheet}" >&2
        exit 1
    fi

    # Verify that the output directory exist, if not create it
    if [ ! -d "${qc_output_dir}" ]; then
        echo "[WARNING] Output directory ${qc_output_dir} does not exist. Creating it."
        mkdir -p "${qc_output_dir}"
    fi

    ##########################################
    #                   Logging              #
    ##########################################

    echo "[INFO] start: \$(date '+%F %T')"
    echo "[INFO] run_id=${run_id}"
    echo "[INFO] bcl_dir=${bcl_dir}"
    echo "[INFO] sample_sheet=${preprocessed_sample_sheet}"
    echo "[INFO] cpus=${task.cpus} mem_gb=${task.memory.toGiga()}"


    ##########################################
    #                   Run                  #
    ##########################################

    # See if run_ID is FlowCellID HCHNTDMX2 or CelineID 260323_A01789_0447_AHCHNTDMX2 (I think FlowCellID HCHNTDMX2 is better)
    cellranger mkfastq \\
        --id="${run_id}" \\ 
        --run="${bcl_dir}" \\
        --csv="${preprocessed_sample_sheet}" \\
        --output-dir=fastq_output \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} 

    EXIT_CODE=$?

    echo "[INFO] end: \$(date '+%F %T')"

    # -----------------------------------------------------------------------
    # Verify successful execution
    # -----------------------------------------------------------------------
    if [ \$EXIT_CODE -ne 0 ]; then
        echo "[ERROR] cellranger mkfastq failed with exit code \$EXIT_CODE." | tee -a logs/mkfastq_${run_id}.log
        exit \$EXIT_CODE
    fi

    # Ensure at least one FASTQ file was generated
    FASTQ_COUNT=\$(find fastq_output -name "*.fastq.gz" | wc -l)
    if [ "\$FASTQ_COUNT" -eq 0 ]; then
        echo "[ERROR] No FASTQ files generated in fastq_output/" | tee -a logs/mkfastq_${run_id}.log
        exit 1
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] mkfastq completed - \$FASTQ_COUNT FASTQ files generated." \\
        | tee -a logs/mkfastq_${run_id}.log

    # -----------------------------------------------------------------------
    # Record tool version
    # -----------------------------------------------------------------------
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
    END_VERSIONS
    """
}
