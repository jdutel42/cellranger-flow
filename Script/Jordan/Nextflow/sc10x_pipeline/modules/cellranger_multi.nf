/*
========================================================================================
    MODULE : CELLRANGER_MULTI
========================================================================================
    Description : Runs one Cell Ranger Multi task per batch.
                  With executor=slurm, each task is submitted as one SLURM job,
                  enabling parallel alignment across all requested batches.
    Outil       : Cell Ranger 7.2.0
    Conda       : /labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_
----------------------------------------------------------------------------------------
    Inputs :
        tuple val(run_id), val(batch)                       → run identifier + batch name
        path(ref_gex)                                       → GEX reference
        path(ref_vdj)                                       → VDJ reference
        val(path_fastq)                                     → FASTQ base directory
        val(fastq_folder_gex)                               → GEX FASTQ subfolder
        val(fastq_folder_vdj)                               → VDJ FASTQ subfolder
        val alignment_output_dir                            → alignment output directory
        val alignment_log_dir                               → alignment log directory
    Outputs :
        tuple val(batch), path("multi_output/${batch}/outs/filtered_feature_bc_matrix/") → ch_matrices
        tuple val(batch), path("multi_output/${batch}/outs/metrics_summary.csv") → ch_metrics
        tuple val(batch), path("multi_output/${batch}/outs/web_summary.html") → ch_web_summaries
        tuple val(batch), path("multi_output/${batch}/outs/molecule_info.h5") → ch_molecule_info
        path "versions.yml" → ch_versions
        path "logs/*.log" → ch_logs
========================================================================================
*/

process CELLRANGER_MULTI {

    // Tag and label for logging and resource allocation
    tag "multi | batch: ${batch}"
    label 'process_high'

    // Use conda environment for reproducibility
    conda '/labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_'

    // Publish key outputs to output_dir
    publishDir (
        path    : alignment_output_dir, // Use the alignment output directory for Cellranger Multi results
        mode    : 'copy',
        pattern : "multi_output/**",
        saveAs  : { filename -> filename }
    )
    publishDir (
        path    : alignment_log_dir, // Use the alignment log directory for Cellranger Multi logs
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    input:
        tuple val(run_id), val(batch)
        path ref_gex
        path ref_vdj
        val path_fastq
        val fastq_folder_gex
        val fastq_folder_vdj
        val alignment_output_dir // Output directory for Cellranger Multi results (val because it's STRING path)
        val alignment_log_dir // Log directory for Cellranger Multi logs (val because it's STRING)

    output:
        tuple val(batch), path("multi_output/${batch}/outs/filtered_feature_bc_matrix/"), emit: matrices
        tuple val(batch), path("multi_output/${batch}/outs/metrics_summary.csv"),         emit: metrics
        tuple val(batch), path("multi_output/${batch}/outs/web_summary.html"),            emit: web_summaries
        tuple val(batch), path("multi_output/${batch}/outs/molecule_info.h5"),            emit: molecule_info
        path "versions.yml",                                                                  emit: versions
        path "logs/*.log",                                                                    emit: logs

    script:
        """
        set -euo pipefail

        mkdir -p logs multi_output

        FASTQ_GEX_DIR="${path_fastq}/${fastq_folder_gex}"
        FASTQ_VDJ_DIR="${path_fastq}/${fastq_folder_vdj}"

        if [ ! -d "${path_fastq}" ]; then
            echo "ERROR: FASTQ base directory does not exist: ${path_fastq}" >&2
            exit 1
        fi

        if [ ! -d "\$FASTQ_GEX_DIR" ]; then
            echo "ERROR: GEX FASTQ directory does not exist: \$FASTQ_GEX_DIR" >&2
            exit 1
        fi

        if [ ! -d "\$FASTQ_VDJ_DIR" ]; then
            echo "ERROR: VDJ FASTQ directory does not exist: \$FASTQ_VDJ_DIR" >&2
            exit 1
        fi

        if ! find "\$FASTQ_GEX_DIR" -maxdepth 1 -name "${batch}_GEX*.fastq.gz" | grep -q .; then
            echo "ERROR: Missing GEX FASTQ files for ${batch} in \$FASTQ_GEX_DIR" >&2
            exit 1
        fi

        if ! find "\$FASTQ_VDJ_DIR" -maxdepth 1 -name "${batch}_VDJ*.fastq.gz" | grep -q .; then
            echo "ERROR: Missing VDJ FASTQ files for ${batch} in \$FASTQ_VDJ_DIR" >&2
            exit 1
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start cellranger multi for batch ${batch}" \
            | tee logs/multi_${batch}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_id=${run_id}" | tee -a logs/multi_${batch}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] cpus=${task.cpus} mem_gb=${task.memory.toGiga()}" | tee -a logs/multi_${batch}.log

        cat > "config_sample_${batch}.csv" <<EOF
[gene-expression]
ref,${ref_gex}
no-bam,FALSE
no-secondary,FALSE

[vdj]
ref,${ref_vdj}

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
${batch}_GEX,\$FASTQ_GEX_DIR,any,${batch}_GEX,Gene Expression,
${batch}_VDJ,\$FASTQ_VDJ_DIR,any,${batch}_VDJ,VDJ,
EOF

        cellranger multi \
            --id="${batch}" \
            --csv="config_sample_${batch}.csv" \
            --localcores=${task.cpus} \
            --localmem=${task.memory.toGiga()} \
            2>&1 | tee -a logs/multi_${batch}.log

        EXIT_CODE=\${PIPESTATUS[0]}

        if [ \$EXIT_CODE -ne 0 ]; then
            echo "[ERROR] cellranger multi failed (code \$EXIT_CODE) for ${batch}." | tee -a logs/multi_${batch}.log
            exit \$EXIT_CODE
        fi

        mv "${batch}" "multi_output/${batch}"

        REQUIRED_OUTPUTS=(
            "multi_output/${batch}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
            "multi_output/${batch}/outs/metrics_summary.csv"
            "multi_output/${batch}/outs/web_summary.html"
            "multi_output/${batch}/outs/molecule_info.h5"
        )

        for output_file in "\${REQUIRED_OUTPUTS[@]}"; do
            if [ ! -f "\$output_file" ] && [ ! -d "\$output_file" ]; then
                echo "[ERROR] Missing expected output: \$output_file" | tee -a logs/multi_${batch}.log
                exit 1
            fi
        done

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] cellranger multi completed for ${batch}." \
            | tee -a logs/multi_${batch}.log

        # -----------------------------------------------------------------------
        # Record tool version
        # -----------------------------------------------------------------------
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
        END_VERSIONS
        """
}
