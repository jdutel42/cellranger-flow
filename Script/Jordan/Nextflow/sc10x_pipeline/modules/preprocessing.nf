/*
========================================================================================
    MODULE : PREPROCESSING
========================================================================================
    Description : Take a raw sample sheet (CSV) and perform necessary preprocessing steps to
                  prepare it for downstream analysis. This may include:
                    - Validating required columns (e.g., sample_id, lane, index)
                    - Normalizing column names and formats
                    - Generating a standardized sample sheet for Cell Ranger
----------------------------------------------------------------------------------------
    Inputs :
        path(raw_sample_sheet)       -> Raw sample sheet (CSV)
    Outputs :
        path "Index_mkfastq_${run_id}.csv"   → ch_standardized_sheet
        path "1_preprocessing_versions.yml"   → ch_versions
========================================================================================
*/

process PREPROCESSING {

    // Tag and label for logging and resource allocation
    // A task is 1 execution of 1 process for 1 set of inputs. for exemple, cellranger_multi with run_id = HCHNTDMX2 and batch_id = 74 is a task
    // A tag is useful to differentiate tasks of the same process (e.g. cellranger_multi) with different inputs (e.g. batch_id) and thus to have different log files and job names in SLURM. The tag is defined in the modules.
    tag "${task.process.toLowerCase()}_${run_id}" // Log the name of the input sample sheet
    label 'process_high' // Use a high resource label since this is a lightweight step

    // Publish the standardized sample sheet to the output directory for reference
    // Publish the standardized sample sheet with a descriptive name
    publishDir (
        path    : "${params.preprocessing_output_dir}", // Use the output directory defined in params 
        mode    : 'copy',
        pattern : "preprocessed_sample_sheet.csv",
        saveAs  : { _filename -> "Index_mkfastq_${run_id}.csv" }
    )
    // Publish logs to a dedicated directory
    publishDir (
        path    : "${params.preprocessing_log_dir}",
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    // Declare process inputs
    input:
    path raw_sample_sheet_file_path // raw_sample_sheet.csv path define in the main.nf file defined in the command line
    val run_id // run_id from params to use in logging and output naming
    val preprocessing_output_dir // Output directory for standardized sample sheet
    val preprocessing_log_dir // Log directory for preprocessing logs
    val today_date // Today's date for logging and output naming

    // Output the standardized sample sheet and versions information
    output:
    path "Index_mkfastq_${run_id}.csv", emit: preprocessed_sample_sheet // Name of the standardized sample sheet channel
    path "1_preprocessing_versions.yml", emit: versions // Emit versions information for reproducibility

    // Script section to perform the preprocessing
    script:
    """
    # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)
    set -euo pipefail

    # Create logs directory if it doesn't exist
    mkdir -p logs

    # Redirect all log (stdout and stderr) to a log file for this process
    exec > >(tee logs/${today_date}_preprocessing_${run_id}.log) 2>&1

    # Load shared logging helpers
    source "${params.logging_script}"

    log_info "
    ╔═══════════════════════════════════════════════════════════════════════════════╗
    ║                         Cell Ranger mkfastq Process Script                    ║
    ╠═══════════════════════════════════════════════════════════════════════════════╣
    ║ Logging input parameters:
    ║ - run_id: ${run_id}
    ║ - raw_sample_sheet: ${raw_sample_sheet_file_path}
    ║ - cpus: ${params.cpu_limit}
    ║ - mem_gb: ${params.memory_limit}
    ║ - preprocessing_output_dir: ${preprocessing_output_dir}
    ║ - preprocessing_log_dir: ${preprocessing_log_dir}
    ║ - today_date: ${today_date}
    ╚═══════════════════════════════════════════════════════════════════════════════╝
    "

    log_start "Starting sample sheet preprocessing: ${raw_sample_sheet_file_path}..."

    # ============================================================================
    # Initialisation
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
    mv "\${TEMP_SHEET}.tmp" Index_mkfastq_${run_id}.csv

    log_ok "Processed columns and cleaned up sample sheet."
    
    log_save "Output saved to Index_mkfastq_${run_id}.csv"

    # -----------------------------------------------------------------------
    # Record tool version
    # -----------------------------------------------------------------------
    log_start "Recording tool versions for reproducibility..."

    cat <<EOF > 1_preprocessing_versions.yml
    "${task.process}":
    awk: "\$(awk -W version 2>&1 | head -n1 | sed 's/\\t/ /g')"
    EOF
        log_ok "Tool versions recorded successfully in 1_preprocessing_versions.yml"

    # -----------------------------------------------------------------------
    # End
    # -----------------------------------------------------------------------
    
    log_save "Preprocessed sample sheet for run_id ${run_id} saved to ${preprocessing_output_dir}/Index_mkfastq_${run_id}.csv."
    log_log "Versions information will be saved to ${params.log_dir}/versions.yml"
    log_log "Logs saved to ${preprocessing_log_dir}/${today_date}_preprocessing_${run_id}.log"

    log_success "Preprocessing of sample sheet completed successfully : ${preprocessing_output_dir}/Index_mkfastq_${run_id}.csv !"
    """

    stub:
        """
        mkdir -p logs
        
        # Create minimal CSV stub with headers
        echo "Sample,Index" > Index_mkfastq_${run_id}.csv
        echo "test_sample,AAAAAA" >> Index_mkfastq_${run_id}.csv
        
        # Create minimal versions file
        cat > 1_preprocessing_versions.yml <<EOF
"PREPROCESSING":
        "run_id": "${run_id}"
        "timestamp": "${today_date}"
EOF
    """
}
