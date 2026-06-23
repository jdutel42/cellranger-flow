/*
========================================================================================
    MODULE : PROCESSING_SAMPLE_SHEET
========================================================================================
    Description : Take a raw sample sheet (CSV) and perform necessary processing steps to
                  prepare it for downstream analysis. This may include:
                    - Validating required columns (e.g., sample_id, lane, index)
                    - Normalizing column names and formats
                    - Generating a processed sample sheet for Cell Ranger
----------------------------------------------------------------------------------------
    Inputs :
        path(raw_sample_sheet_file_path) -> Raw sample sheet path for input

    Outputs :
        path "Index_mkfastq_${params.bcl_id}.csv"         → ch_processed_sample_sheet
        path "Processing_versions_${params.bcl_id}.yml"   → ch_versions
========================================================================================
*/

process PROCESSING_SAMPLE_SHEET {

    // Tag and label for logging and resource allocation
    // A task is 1 execution of 1 process for 1 set of inputs. for exemple, cellranger_multi with run_id = HCHNTDMX2 and batch_id = 74 is a task
    // A tag is useful to differentiate tasks of the same process (e.g. cellranger_multi) with different inputs (e.g. batch_id) and thus to have different log files and job names in SLURM. The tag is defined in the modules.
    tag "PrSS_${params.bcl_id}_Processing_Sample_Sheet_${params.protocol_id}" // Log the name of the input sample sheet for "Pro" = "Processing sample-sheet"
    label 'process_high' // Use a high resource label since this is a lightweight step

    // Publish the processed sample sheet to the output directory for reference
    publishDir (
        path    : "${params.processing_output_dir}", // Use the output directory defined in params 
        mode    : 'copy',
        pattern : "Index_mkfastq_${params.bcl_id}.csv", // Publish the processed sample sheet with a descriptive name
        saveAs  : { _filename -> "Index_mkfastq_${params.bcl_id}.csv" }
    )

    // Additionally, publish the processed sample sheet to a specific directory for traceability
    publishDir (
        path    : "/home/dutel/Data/Sample_sheet/Modified",
        mode    : 'copy',
        pattern : "Index_mkfastq_${params.bcl_id}.csv",
        saveAs  : { _filename -> "Index_mkfastq_${params.bcl_id}.csv" }
    )

    // Publish logs to a dedicated directory
    publishDir (
        path    : "${params.processing_log_dir}",
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    // Declare process inputs
    input:
    path raw_sample_sheet_file_path // raw_sample_sheet.csv path define in the main.nf file defined in the command line

    // Output the processed sample sheet and versions information
    output:
    path "Index_mkfastq_${params.bcl_id}.csv",       emit: ch_processed_sample_sheet // Name of the processed sample sheet channel
    path "Processing_versions_${params.bcl_id}.yml", emit: ch_versions // Emit versions information for reproducibility

    // Script section to perform the processing
    script:
    """
    # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)
    set -euo pipefail

    # Create processing_output and logs directories in work/
    mkdir -p processing_output_${params.bcl_id} logs

    # Redirect all logs (stdout and stderr) to a log file for this process
    exec > >(tee logs/Processing_${params.bcl_id}.log) 2>&1

    # Load shared logging helpers
    source "${params.logging_script}"

    log_init "Step 1: Processing sample sheet = ${raw_sample_sheet_file_path}..."

    log_log "Logs will be saved to ${params.processing_log_dir}/Processing_${params.bcl_id}.log"

    log_info "
    ╔═══════════════════════════════════════════════════════════════════════════════╗
    ║                         Processing Process Script                             ║
    ╠═══════════════════════════════════════════════════════════════════════════════╣
    ║ Logging input parameters:
    ║ - bcl_id: ${params.bcl_id}
    ║ - run_id: ${params.run_id}
    ║ - raw_sample_sheet: ${params.raw_sample_sheet_file_path}
    ║ - cpus: ${params.cpu_limit}
    ║ - mem_gb: ${params.memory_limit}
    ║ - processing_output_dir: ${params.processing_output_dir}
    ║ - processing_log_dir: ${params.processing_log_dir}
    ╚═══════════════════════════════════════════════════════════════════════════════╝
    "

    # ============================================================================
    # Initialisation & Verification
    # ============================================================================

    # Create temporary file for processing
    TEMP_SHEET=\$(mktemp)
    trap "rm -f \$TEMP_SHEET \${TEMP_SHEET}.tmp" EXIT

    # ============================================================================
    # Run awk
    # ============================================================================
    log_start "Removing header, reads, settings sections and [Data] title..."

    # Remove [Header], [Reads], and [Settings] sections
    # Keep only [Data] section and remove the [Data] header
    awk '
        BEGIN { in_data = 0 }
        /^\\[Data\\]/ { in_data = 1; next }
        /^\\[Header\\]/ || /^\\[Reads\\]/ || /^\\[Settings\\]/ { in_data = 0; next }
        in_data && NF > 0 { print }
        ' "${raw_sample_sheet_file_path}" > "\$TEMP_SHEET"

    log_ok "Extracted [Data] section from raw sample sheet"

    log_start "Processing columns..."

    # Process columns
    # 1. Replace Species with Lane (*)
    # 2. Remove Sample_Name
    # 3. Rename Sample_ID → Sample
    # 4. Rename index10X → Index
    # 5. Remove quotes

    awk -F',' -v OFS=',' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
            col = \$i
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", col)
            gsub(/"/, "", col)

            if (col == "Sample_ID")      header[i] = "Sample"
            else if (col == "index10X")  header[i] = "Index"
            else if (col == "Species") { lane_idx = i; header[i] = "Lane" }
            else if (col == "Sample_Name") skip[i] = 1
            else header[i] = col
            }

            out = ""
            for (i = 1; i <= NF; i++) {
            if (!(i in skip)) {
                if (out != "") out = out OFS
                out = out header[i]
            }
            }
            print out
            next
        }

        {
            out = ""
            for (i = 1; i <= NF; i++) {
            if (i in skip) continue

            col = \$i
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", col)
            gsub(/"/, "", col)

            if (i == lane_idx) col = "*"

            if (out != "") out = out OFS
            out = out col
            }
            print out
        }
    ' "\$TEMP_SHEET" > "\${TEMP_SHEET}.tmp"

    # ============================================================================
    # Verification & End
    # ============================================================================
    log_verify "Verifying processed sample sheet..."

    if [ ! -s "\${TEMP_SHEET}.tmp" ]; then
        log_error "Processed sample sheet is empty. Check the input file and preprocessing steps."
    fi

    # Rename the processed sheet to the final output name
    mv "\${TEMP_SHEET}.tmp" Index_mkfastq_${params.bcl_id}.csv

    log_ok "Processed columns and cleaned up sample sheet."
    
    log_save "Output saved to Index_mkfastq_${params.bcl_id}.csv"

    # -----------------------------------------------------------------------
    # Record tool version
    # -----------------------------------------------------------------------
    log_start "Recording tool versions for reproducibility..."

    cat <<EOF > Processing_versions_${params.bcl_id}.yml
    "${task.process}":
    awk: "\$(awk -W version 2>&1 | head -n1 | sed 's/\\t/ /g')"
    EOF
    
    log_ok "Tool versions recorded successfully in Processing_versions_${params.bcl_id}.yml"

    # -----------------------------------------------------------------------
    # End
    # -----------------------------------------------------------------------
    
    log_save "Processed sample sheet for bcl_id ${params.bcl_id} saved to ${params.processing_output_dir}/Index_mkfastq_${params.bcl_id}.csv."
    log_log "Versions information will be saved to ${params.run_traceability_log_dir}/Processing_versions_${params.bcl_id}.yaml"
    log_log "Logs saved to ${params.processing_log_dir}/Processing_${params.bcl_id}.log"

    log_success "Processing of sample sheet completed successfully : ${params.processing_output_dir}/Index_mkfastq_${params.bcl_id}.csv !"
    """

    stub:
        """
        mkdir -p logs
        
        # Create minimal CSV stub with headers
        echo "Sample,Index" > Index_mkfastq_${params.bcl_id}.csv
        echo "test_sample,AAAAAA" >> Index_mkfastq_${params.bcl_id}.csv
        
        # Create minimal versions file
        cat > Processing_versions_${params.bcl_id}.yml <<EOF
"PROCESSING":
        "bcl_id": "${params.bcl_id}"
        "run_id": "${params.run_id}"
        "awk_version": "awk 5.1.0"
EOF
    """
}
