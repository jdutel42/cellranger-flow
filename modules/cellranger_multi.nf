/*
========================================================================================
    MODULE : CELLRANGER_MULTI
========================================================================================
    Description : Runs one Cell Ranger Multi task per sample.
                  With executor=slurm, each task is submitted as one SLURM job,
                  enabling parallel alignment across all requested batches.
    Outil       : Cell Ranger 7.2.0
    Conda       : /labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_
----------------------------------------------------------------------------------------
    Inputs :
        tuple path(fastq_dir), val(sample_id)

    Outputs :
        tuple val(sample_id), path("multi_output/${sample_id}/outs/metrics_summary.csv") -> ch_metrics
        tuple val(sample_id), path("multi_output/${sample_id}/outs/web_summary.html") -> ch_web_summaries
        path "Cellranger_multi_versions_${params.bcl_id}_${sample_id}.yml" -> ch_versions
========================================================================================
*/

process CELLRANGER_MULTI {

    // Tag and label for logging and resource allocation
    tag "CrMu_${sample_id}_${params.bcl_id}_CellRanger_Multi_${params.protocol_prefix}" // Log the run ID for traceability for "CrM" = "Cell Ranger Multi"
    label 'process_high'

    // Use conda environment for reproducibility
    conda '/labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_'

    // Publish the full batch output tree to the alignment results directory
    publishDir (
        path    : params.alignment_output_dir, // Use the alignment output directory for Cellranger Multi results
        mode    : 'copy',
        pattern : "multi_output_${params.bcl_id}/${sample_id}/**"
    )
    
    publishDir (
        path    : params.alignment_log_dir, // Use the alignment log directory for Cellranger Multi logs
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    input:
        tuple path(fastq_dir), val(sample_id) // Channel of tuples containing (fastq_dir, sample_id) for each sample

    output:
        tuple val(sample_id), path("multi_output_${params.bcl_id}/${sample_id}"),                                                               emit: ch_multi_output
        tuple val(sample_id), path("multi_output_${params.bcl_id}/${sample_id}/outs/per_sample_outs/${sample_id}/metrics_summary.csv"),         emit: ch_metrics
        tuple val(sample_id), path("multi_output_${params.bcl_id}/${sample_id}/outs/per_sample_outs/${sample_id}/web_summary.html"),            emit: ch_web_summaries
        path "Cellranger_multi_versions_${params.bcl_id}_${sample_id}.yml",                                                                     emit: ch_versions

    script:
        """
        # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)
        set -euo pipefail

        # Create multi_output and logs directories in work/
        mkdir -p logs multi_output_${params.bcl_id}/${sample_id}

        # Redirect all log (stdout and stderr) to a log file for this process
        exec > >(tee -a logs/Cellranger_multi_log_${params.bcl_id}_${sample_id}.log) 2>&1

        # Load shared logging helpers
        source "${params.logging_script}"

        log_init "Step 3: Performing alignment with Cellranger Multi for run_id = ${params.run_id} & sample_id = ${sample_id}..."

        log_log "Logs will be saved to ${params.alignment_log_dir}/Cellranger_multi_log_${params.bcl_id}_${sample_id}.log"

        log_info "
        ╔═══════════════════════════════════════════════════════════════════════════════╗
        ║                         Cell Ranger multi Process Script                      ║
        ╠═══════════════════════════════════════════════════════════════════════════════╣
        ║ Logging input parameters:
        ║ - run_id: ${params.run_id}
        ║ - bcl_id: ${params.bcl_id}
        ║ - sample_id: ${sample_id}
        ║ - fastq_folder_gex: ${params.fastq_folder_gex}
        ║ - fastq_folder_vdj: ${params.fastq_folder_vdj}
        ║ - fastq_folder_adt: ${params.fastq_folder_adt}
        ║ - genome_reference_path: ${params.genome_reference_path}
        ║ - vdj_reference_path: ${params.vdj_reference_path}
        ║ - cpus: ${params.cpu_limit}
        ║ - mem_gb: ${params.memory_limit}
        ║ - alignment_output_dir: ${params.alignment_output_dir}
        ║ - alignment_log_dir: ${params.alignment_log_dir}
        ╚═══════════════════════════════════════════════════════════════════════════════╝
        "

        # -------------------------------------------------------------------------
        # Initialisation & Verification
        # -------------------------------------------------------------------------
        
        log_verify "Verifying input files and directories..."

        # Dynamically set FASTQ directories based on provided parameters or default to the fastq_dir from the channel
        # If fastq_folder_gex, fastq_folder_vdj, or fastq_folder_adt are not provided, default fastq_dir to the provided fastq_dir from the channel
        if [ -z "${params.fastq_folder_gex}" ]; then
            FASTQ_GEX_DIR="${fastq_dir}"
        else
            FASTQ_GEX_DIR="${params.fastq_folder_gex}"
        fi
        if [ -z "${params.fastq_folder_vdj}" ]; then
            FASTQ_VDJ_DIR="${fastq_dir}"
        else
            FASTQ_VDJ_DIR="${params.fastq_folder_vdj}"
        fi
        if [ -z "${params.fastq_folder_adt}" ]; then
            FASTQ_ADT_DIR="${fastq_dir}"
        else
            FASTQ_ADT_DIR="${params.fastq_folder_adt}"
        fi

        # Verification of input files and directories
        if [ ! -d "\${FASTQ_GEX_DIR}" ]; then
            log_error "FASTQ GEX directory does not exist: \${FASTQ_GEX_DIR}"
        fi
        if [ ! -d "\${FASTQ_VDJ_DIR}" ]; then
            log_error "FASTQ VDJ directory does not exist: \${FASTQ_VDJ_DIR}"
        fi
        if [ ! -d "\${FASTQ_ADT_DIR}" ]; then
            log_error "FASTQ ADT directory does not exist: \${FASTQ_ADT_DIR}"
        fi

        # Verify that the FASTQ directories contain at least one FASTQ file
        # Note: Uncomment the following checks if you want to enforce the presence of FASTQ files
        #if [[ -z \$(find "\$FASTQ_GEX_DIR" -type f -name "*.fastq.gz" -print -quit) ]]; then
        #    log_error "Missing GEX FASTQ files for ${sample_id} in \$FASTQ_GEX_DIR"
        #fi
        #if [[ -z \$(find "\$FASTQ_VDJ_DIR" -type f -name "*.fastq.gz" -print -quit) ]]; then
        #    log_error "Missing VDJ FASTQ files for ${sample_id} in \$FASTQ_VDJ_DIR"
        #fi
        #if [[ -z \$(find "\$FASTQ_ADT_DIR" -type f -name "*.fastq.gz" -print -quit) ]]; then
        #    log_error "Missing ADT FASTQ files for ${sample_id} in \$FASTQ_ADT_DIR"
        #fi

        log_ok "Input files and directories verified successfully."

        # -------------------------------------------------------------------------
        # Create Cell Ranger multi config CSV for this batch
        # -------------------------------------------------------------------------
        
        log_start "Creating Cell Ranger multi config CSV for run_id ${params.run_id} & sample_id ${sample_id}..."

        # Build the config file incrementally depending on requested libraries
        > "config_sample_${sample_id}.csv"

        if [ ${params.gex} ]; then
            cat >> "config_sample_${sample_id}.csv" <<EOF
[gene-expression]
ref,${params.genome_reference_path}
no-bam,FALSE
no-secondary,FALSE
EOF
        fi

        if [ ${params.vdj} ]; then
            cat >> "config_sample_${sample_id}.csv" <<EOF
[vdj]
ref,${params.vdj_reference_path}
EOF
        fi

        if [ ${params.adt} ]; then
            cat >> "config_sample_${sample_id}.csv" <<EOF
[feature]
ref,${params.adt_reference_path}
EOF
        fi

        # Libraries section: always present, but entries depend on enabled libraries
        cat >> "config_sample_${sample_id}.csv" <<EOF
[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
EOF

        if [ ${params.gex} ]; then
            echo "${sample_id}_GEX,\$FASTQ_GEX_DIR,any,${sample_id}_GEX,Gene Expression," >> "config_sample_${sample_id}.csv"
        fi

        if [ ${params.vdj} ]; then
            echo "${sample_id}_VDJ,\$FASTQ_VDJ_DIR,any,${sample_id}_VDJ,VDJ," >> "config_sample_${sample_id}.csv"
        fi

        if [ ${params.adt} ]; then
            echo "${sample_id}_ADT,\$FASTQ_ADT_DIR,any,${sample_id}_ADT,Antibody Capture," >> "config_sample_${sample_id}.csv"
        fi

        # Add [samples] section with hashtag mappings if ADT is enabled
        if [ ${params.adt} ] && [ -n "${params.adt_samples_hashtags}" ] && [ "${params.adt_samples_hashtags}" != "{}" ]; then
            cat >> "config_sample_${sample_id}.csv" <<EOF
[samples]
sample_id,hashtag_ids
${params.adt_samples_hashtags.collect {label, hashtag_id -> "${label},${hashtag_id}" }.join('\n')}
EOF
        fi

        log_ok "Config CSV created: config_sample_${sample_id}.csv"

        log_save "Config CSV content saved in logs/config_sample_${sample_id}.csv. Add a publishDir if you want to save it to output directories."

        # -------------------------------------------------------------------------
        # Run Cellranger multi
        # -------------------------------------------------------------------------
        
        log_start "Running Cell Ranger multi for run_id ${params.run_id} & sample_id ${sample_id}..."
        
        /labos/UGM/dev/${params.cellranger_version}/bin/cellranger multi \
            --id="${sample_id}" \
            --csv="config_sample_${sample_id}.csv" \
            --localcores=${params.cpu_limit} \
            --localmem=${params.memory_limit} \
            2>&1 | tee -a logs/Cellranger_multi_log_${params.bcl_id}_${sample_id}.log

        log_ok "Cell Ranger multi finished for run_id ${params.run_id} & sample_id ${sample_id}."

        # -----------------------------------------------------------------------
        # Verify successful execution
        # -----------------------------------------------------------------------

        log_verify "Verifying cellranger multi output for ${sample_id}..."

        mv "${sample_id}" "multi_output/${sample_id}"

        # Verify key outputs exist
        if [ ! -f "multi_output/${sample_id}/outs/per_sample_outs/${sample_id}/metrics_summary.csv" ]; then
            log_error "Missing metrics_summary.csv"
        fi

        if [ ! -f "multi_output/${sample_id}/outs/per_sample_outs/${sample_id}/web_summary.html" ]; then
            log_error "Missing web_summary.html"
        fi

        log_ok "Cellranger multi completed successfully for ${sample_id}."

        log_save "Cellranger multi outputs (metrics_summary.csv, web_summary.html, filtered_feature_bc_matrix...) saved in ${params.alignment_output_dir}."

        # -----------------------------------------------------------------------
        # Record tool version
        # -----------------------------------------------------------------------
        
        log_start "Recording tool versions for reproducibility..."
        
        cat <<-END_VERSIONS > Cellranger_multi_versions_${params.bcl_id}_${sample_id}.yml
        "${task.process}":
            cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
        END_VERSIONS

        log_ok "Tool versions recorded successfully in Cellranger_multi_versions_${params.bcl_id}_${sample_id}.yml"

        # -----------------------------------------------------------------------
        # End
        # -----------------------------------------------------------------------

        log_save "Cellranger multi output fo run_id ${params.run_id} & sample_id ${sample_id} saved to ${params.alignment_output_dir}/multi_output/${sample_id}/."
        log_log "Versions information will be saved to ${params.run_traceability_log_dir}/Cellranger_multi_versions_${params.bcl_id}_${sample_id}.yml"

        log_log "Logs saved to ${params.alignment_log_dir}/Cellranger_multi_log_${params.bcl_id}_${sample_id}.log"

        log_success "CellRanger multi process completed successfully for run_id = ${params.run_id} & sample_id = ${sample_id} and results are available at ${params.alignment_output_dir}!"
        """

    stub:
    """
    mkdir -p multi_output/${sample_id}/outs logs
    
    # Create minimal metrics CSV
    cat > multi_output/${sample_id}/outs/metrics_summary.csv <<EOF
Metric,Value
"Estimated Number of Cells","1000"
"Mean Reads per Cell","50000"
"Median UMI Counts per Cell","5000"
EOF
    
    # Create minimal HTML report
    cat > multi_output/${sample_id}/outs/web_summary.html <<EOF
<!DOCTYPE html>
<html>
<head><title>Cell Ranger Multi Report - ${sample_id}</title></head>
<body>
<h1>Cell Ranger Multi Report (STUB)</h1>
<p>Run ID: ${params.run_id}</p>
<p>Sample ID: ${sample_id}</p>
<p>Estimated Cells: 1000</p>
</body>
</html>
EOF
    
    # Create minimal versions file
    cat > Cellranger_multi_versions_${params.bcl_id}_${sample_id}.yml <<EOF
"CELLRANGER_MULTI":
    "cellranger": "${params.cellranger_version}"
    "sample_id": "${sample_id}"
EOF
    """
}
