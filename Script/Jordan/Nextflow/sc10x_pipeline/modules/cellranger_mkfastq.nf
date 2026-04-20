/*
========================================================================================
    MODULE : CELLRANGER_MKFASTQ
========================================================================================
    Description : Converts a BCL directory into FASTQ files using cellranger mkfastq.
    Outil       : Cell Ranger 7.2.0
    Container   : nfcore/cellranger:7.2.0
----------------------------------------------------------------------------------------
    Inputs :
        tuple val(run_id), path(bcl_dir), path(sample_sheet)
    Outputs :
        tuple val(run_id), path("fastq_output/*")   → ch_fastqs
        path "versions.yml"                          → ch_versions
        path "logs/*.log"                            → ch_logs
========================================================================================
*/

process CELLRANGER_MKFASTQ {

    tag "mkfastq | run: ${run_id}"
    label 'process_high'

    container 'nfcore/cellranger:7.2.0'

    publishDir (
        path    : "${params.output_dir}/mkfastq/${run_id}",
        mode    : 'copy',
        pattern : "fastq_output/**",
        saveAs  : { filename -> filename }
    )
    
    publishDir (
        path    : "${params.output_dir}/logs/mkfastq",
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    input:
        tuple val(run_id), path(bcl_dir), path(sample_sheet)

    output:
        tuple val(run_id), path("fastq_output/*"), emit: fastqs
        path "versions.yml",                       emit: versions
        path "logs/*.log",                         emit: logs

    when:
        task.ext.when == null || task.ext.when

    script:
        // Build optional arguments
        def args = task.ext.args ?: ''

        """
        # -----------------------------------------------------------------------
        # Input validation
        # -----------------------------------------------------------------------
        if [ ! -d "${bcl_dir}" ]; then
            echo "ERROR: BCL directory does not exist: ${bcl_dir}" >&2
            exit 1
        fi

        if [ ! -f "${sample_sheet}" ]; then
            echo "ERROR: Sample sheet does not exist: ${sample_sheet}" >&2
            exit 1
        fi

        # -----------------------------------------------------------------------
        # Create output and log directories
        # -----------------------------------------------------------------------
        mkdir -p fastq_output logs

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting cellranger mkfastq - run_id: ${run_id}" \\
            | tee logs/mkfastq_${run_id}.log

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] BCL dir    : ${bcl_dir}"   | tee -a logs/mkfastq_${run_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample sheet: ${sample_sheet}" | tee -a logs/mkfastq_${run_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Allocated CPUs: ${task.cpus}" | tee -a logs/mkfastq_${run_id}.log

        # -----------------------------------------------------------------------
        # Run cellranger mkfastq
        # -----------------------------------------------------------------------
        cellranger mkfastq \\
            --id="${run_id}_mkfastq" \\
            --run="${bcl_dir}" \\
            --csv="${sample_sheet}" \\
            --output-dir=fastq_output \\
            --localcores=${task.cpus} \\
            --localmem=${task.memory.toGiga()} \\
            --delete-undetermined \\
            --ignore-dual-index-flowcells \\
            ${args} \\
            2>&1 | tee -a logs/mkfastq_${run_id}.log

        EXIT_CODE=\${PIPESTATUS[0]}

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

    stub:
        """
        mkdir -p fastq_output/Sample_A fastq_output/Sample_B logs

        # Dummy FASTQ files for stub mode
        touch fastq_output/Sample_A/Sample_A_S1_L001_R1_001.fastq.gz
        touch fastq_output/Sample_A/Sample_A_S1_L001_R2_001.fastq.gz
        touch fastq_output/Sample_B/Sample_B_S2_L001_R1_001.fastq.gz
        touch fastq_output/Sample_B/Sample_B_S2_L001_R2_001.fastq.gz
        touch logs/mkfastq_stub.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: "7.2.0"
        END_VERSIONS
        """
}
