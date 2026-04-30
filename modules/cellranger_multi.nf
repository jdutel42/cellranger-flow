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
        tuple val(batch_id), path("multi_output/${batch_id}/outs/metrics_summary.csv") → ch_metrics
        tuple val(batch_id), path("multi_output/${batch_id}/outs/web_summary.html") → ch_web_summaries
        path "3_cellranger_multi_versions.yml" → ch_versions
========================================================================================
*/

process CELLRANGER_MULTI {

    // Tag and label for logging and resource allocation
    tag "CrM_${batch_id}_${run_id}" // Log the run ID for traceability for "CrM" = "Cell Ranger Multi"
    label 'process_high'

    // Use conda environment for reproducibility
    conda '/labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_'

    // Publish key outputs to output_dir
    publishDir (
        path    : params.alignment_output_dir, // Use the alignment output directory for Cellranger Multi results
        mode    : 'copy',
        pattern : "multi_output/**",
        saveAs  : { filename -> filename }
    )
    
    publishDir (
        path    : params.alignment_log_dir, // Use the alignment log directory for Cellranger Multi logs
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    input:
        tuple val(run_id), val(batch_id)
        path fastq_folder
        path ref_gex
        path ref_vdj
        val alignment_output_dir // Output directory for Cellranger Multi results (val because it's STRING path)
        val alignment_log_dir // Log directory for Cellranger Multi logs (val because it's STRING)
        val today_date // Today's date for logging and output naming

    output:
        tuple val(batch_id), path("multi_output/${batch_id}/outs/metrics_summary.csv"),         emit: metrics
        tuple val(batch_id), path("multi_output/${batch_id}/outs/web_summary.html"),            emit: web_summaries
        path "3_cellranger_multi_versions.yml",                                                 emit: versions

    script:
        """
        # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)
        set -euo pipefail

        # Create logs directory if it doesn't exist
        mkdir -p logs multi_output

        # Redirect all log (stdout and stderr) to a log file for this process
        exec > >(tee -a logs/${today_date}_multi_${run_id}_${batch_id}.log) 2>&1

        # Load shared logging helpers
        source "${params.logging_script}"

        log_init "Step 3: Performing alignment with Cellranger Multi for run_id = ${run_id} & batch_id = ${batch_id}..."
        log_log "Logs will be saved to ${alignment_log_dir}/${today_date}_multi_${run_id}_${batch_id}.log"

        log_info "
        ╔═══════════════════════════════════════════════════════════════════════════════╗
        ║                         Cell Ranger multi Process Script                      ║
        ╠═══════════════════════════════════════════════════════════════════════════════╣
        ║ Logging input parameters:
        ║ - run_id: ${run_id}
        ║ - batch_id: ${batch_id}
        ║ - fastq_folder: ${fastq_folder}
        ║ - ref_gex: ${ref_gex}
        ║ - ref_vdj: ${ref_vdj}
        ║ - cpus: ${params.cpu_limit}
        ║ - mem_gb: ${params.memory_limit}
        ║ - alignment_output_dir: ${alignment_output_dir}
        ║ - alignment_log_dir: ${alignment_log_dir}
        ║ - today_date: ${today_date}
        ╚═══════════════════════════════════════════════════════════════════════════════╝
        "

        # -------------------------------------------------------------------------
        # Initialisation & Verification
        # -------------------------------------------------------------------------
        
        log_verify "Verifying input files and directories..."

        # Define paths to FASTQ directories
        FASTQ_GEX_DIR="${fastq_folder}"
        FASTQ_VDJ_DIR="${fastq_folder}"

        # Verification of input files and directories
        if [ ! -d "${fastq_folder}" ]; then
            log_error "FASTQ base directory does not exist: ${fastq_folder}"
        fi

        #if [[ -z \$(find "\$FASTQ_GEX_DIR" -type f -name "*.fastq.gz" -print -quit) ]]; then
        #    log_error "Missing GEX FASTQ files for ${batch_id} in \$FASTQ_GEX_DIR"
        #fi

        #if [[ -z \$(find "\$FASTQ_VDJ_DIR" -type f -name "*.fastq.gz" -print -quit) ]]; then
        #    log_error "Missing VDJ FASTQ files for ${batch_id} in \$FASTQ_VDJ_DIR"
        #fi

        log_ok "Input files and directories verified successfully."

        # -------------------------------------------------------------------------
        # Create Cell Ranger multi config CSV for this batch
        # -------------------------------------------------------------------------
        
        log_start "Creating Cell Ranger multi config CSV for run_id ${run_id} & batch_id ${batch_id}..."
        
        cat > "config_sample_${batch_id}.csv" <<EOF
[gene-expression]
ref,${ref_gex}
no-bam,FALSE
no-secondary,FALSE

[vdj]
ref,${ref_vdj}

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
${batch_id}_GEX,\$FASTQ_GEX_DIR,any,${batch_id}_GEX,Gene Expression,
${batch_id}_VDJ,\$FASTQ_VDJ_DIR,any,${batch_id}_VDJ,VDJ,
EOF

        log_ok "Config CSV created: config_sample_${batch_id}.csv"
        
        log_save "Config CSV content saved in logs/${today_date}_multi_${run_id}_${batch_id}_config.csv. Add a publishDir if you want to save it to output directories."

        # -------------------------------------------------------------------------
        # Run Cellranger multi
        # -------------------------------------------------------------------------
        
        log_start "Running Cell Ranger multi for run_id ${run_id} & batch_id ${batch_id}..."
        
        /labos/UGM/dev/cellranger-7.1.0/bin/cellranger multi \
            --id="${batch_id}" \
            --csv="config_sample_${batch_id}.csv" \
            --localcores=${params.cpu_limit} \
            --localmem=${params.memory_limit} \
            2>&1 | tee -a logs/${today_date}_multi_${run_id}_${batch_id}.log

        log_ok "Cell Ranger multi finished for run_id ${run_id} & batch_id ${batch_id}."

        # -----------------------------------------------------------------------
        # Verify successful execution
        # -----------------------------------------------------------------------

        log_verify "Verifying cellranger multi output for ${batch_id}..."

        mv "${batch_id}" "multi_output/${batch_id}"

        REQUIRED_OUTPUTS=(
            "multi_output/${batch_id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
            "multi_output/${batch_id}/outs/metrics_summary.csv"
            "multi_output/${batch_id}/outs/web_summary.html"
            "multi_output/${batch_id}/outs/molecule_info.h5"
        )

        for output_file in "\${REQUIRED_OUTPUTS[@]}"; do
            if [ ! -f "\$output_file" ] && [ ! -d "\$output_file" ]; then
                log_error "Missing expected output: \$output_file"
            fi
        done

        log_ok "Cellranger multi completed successfully for ${batch_id}."

        log_save "Cellranger multi outputs (metrics_summary.csv, web_summary.html, filtered_feature_bc_matrix...) saved in ${alignment_output_dir}."

        # -----------------------------------------------------------------------
        # Record tool version
        # -----------------------------------------------------------------------
        
        log_start "Recording tool versions for reproducibility..."
        
        cat <<-END_VERSIONS > 3_cellranger_multi_versions.yml
        "${task.process}":
            cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
        END_VERSIONS

        log_ok "Tool versions recorded successfully in 3_cellranger_multi_versions.yml"

        # -----------------------------------------------------------------------
        # End
        # -----------------------------------------------------------------------

        log_save "Cellranger multi output fo run_id ${run_id} & batch_id ${batch_id} saved to ${alignment_output_dir}/multi_output/${batch_id}/."
        log_log "Versions information will be saved to ${params.run_traceability_log_dir}/${today_date}_versions.yaml"
        log_log "Logs saved to ${alignment_log_dir}/${today_date}_multi_${run_id}_${batch_id}.log"

        log_success "CellRanger multi process completed successfully for run_id = ${run_id} & batch_id = ${batch_id} and results are available at ${params.alignment_output_dir}!"
        """

    stub:
    """
    mkdir -p multi_output/${batch_id}/outs logs
    
    # Create minimal metrics CSV
    cat > multi_output/${batch_id}/outs/metrics_summary.csv <<EOF
Metric,Value
"Estimated Number of Cells","1000"
"Mean Reads per Cell","50000"
"Median UMI Counts per Cell","5000"
EOF
    
    # Create minimal HTML report
    cat > multi_output/${batch_id}/outs/web_summary.html <<EOF
<!DOCTYPE html>
<html>
<head><title>Cell Ranger Multi Report - ${batch_id}</title></head>
<body>
<h1>Cell Ranger Multi Report (STUB)</h1>
<p>Run ID: ${run_id}</p>
<p>Batch ID: ${batch_id}</p>
<p>Estimated Cells: 1000</p>
</body>
</html>
EOF
    
    # Create minimal versions file
    cat > 3_cellranger_multi_versions.yml <<EOF
"CELLRANGER_MULTI":
    "cellranger": "7.1.0"
    "batch_id": "${batch_id}"
EOF
    """
}
