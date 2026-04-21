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
        path "versions.yml"                    → ch_versions
========================================================================================
*/

process PREPROCESSING {

    // Tag and label for logging and resource allocation
    tag "preprocess | run_id: ${run_id}" // Log the name of the input sample sheet
    label 'process_low' // Use a low resource label since this is a lightweight step

    // Publish the standardized sample sheet to the output directory for reference
    // Publish the standardized sample sheet with a descriptive name
    publishDir (
        path    : "${preprocessing_output_dir}", // Use the output directory defined in params 
        mode    : 'copy',
        pattern : "preprocessed_sample_sheet.csv",
        saveAs  : { _filename -> "Index_mkfastq_${run_id}.csv" }
    )
    // Publish logs to a dedicated directory
    publishDir (
        path    : "${preprocessing_log_dir}",
        mode    : 'copy',
        pattern : "Pre_processing.log"
    )

    // Declare process inputs
    input:
    path raw_sample_sheet_file_path // raw_sample_sheet.csv path define in the main.nf file defined in the command line
    val run_id // run_id from params to use in logging and output naming
    val preprocessing_output_dir // Output directory for standardized sample sheet
    val preprocessing_log_dir // Log directory for preprocessing logs

    // Output the standardized sample sheet and versions information
    output:
    path "Index_mkfastq_${run_id}.csv", emit: preprocessed_sample_sheet // Name of the standardized sample sheet channel
    path "versions.yml", emit: versions // Emit versions information for reproducibility

    // Script section to perform the preprocessing
    script:
    """

    set -euo pipefail # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)

    exec > >(tee Pre_processing.log) 2>&1

    # ============================================================================
    # Functions
    # ============================================================================

    log_info() {
      echo "[INFO] \$(date '+%Y-%m-%d %H:%M:%S') - \$1"
    }

    log_error() {
      echo "[ERROR] \$(date '+%Y-%m-%d %H:%M:%S') - \$1" >&2
    }

    # ============================================================================
    # Main Script
    # ============================================================================

    log_info "Starting sample sheet preprocessing: ${raw_sample_sheet_file_path}"

    # Create temporary file for processing
    TEMP_SHEET=\$(mktemp)
    trap "rm -f \$TEMP_SHEET \${TEMP_SHEET}.tmp" EXIT

    # Remove [Header], [Reads], and [Settings] sections
    # Keep only [Data] section and remove the [Data] header
    awk '
        BEGIN { in_data = 0 }
        /^\\[Data\\]/ { in_data = 1; next }
        /^\\[Header\\]/ || /^\\[Reads\\]/ || /^\\[Settings\\]/ { in_data = 0; next }
        in_data && NF > 0 { print }
        ' "${raw_sample_sheet_file_path}" > "\$TEMP_SHEET"

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
        
    mv "\${TEMP_SHEET}.tmp" Index_mkfastq_${run_id}.csv
    
    cat <<EOF > versions.yml
"${task.process}":
  awk: "\$(awk -W version 2>&1 | head -n1 | sed 's/\\t/ /g')"
EOF

    log_info "Preprocessing completed: Index_mkfastq_${run_id}.csv"
    """
}
