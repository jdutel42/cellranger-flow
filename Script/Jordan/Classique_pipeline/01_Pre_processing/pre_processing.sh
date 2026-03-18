#!/bin/bash

###############################################################################
# Pre-processing Script for Single Cell Sample Sheet
# Purpose: Prepare and modify SampleSheet for 10X analysis
###############################################################################

set -euo pipefail # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)

# ============================================================================
# Configuration
# ============================================================================

SAMPLE_SHEET_DIR="../../Sample_sheet/Original/"
BCL_BASE_PATH="/media/CRCT13/60To/SingleCell/BclNovaseq"
SHARE_BASE_PATH="/media/CRCT13/20To"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

check_file_exists() {
    if [ ! -f "$1" ]; then
        log_error "File not found: $1"
        exit 1
    fi
}

# ============================================================================
# Main Script
# ============================================================================

main() {
    log_info "Starting Sample Sheet preprocessing..."

    # Check if Sample Sheet directory exists
    if [ ! -d "$SAMPLE_SHEET_DIR" ]; then
        log_error "Sample Sheet directory not found: $SAMPLE_SHEET_DIR"
        exit 1
    fi

    # Get the Sample Sheet file (assuming there's one CSV file)
    SAMPLE_SHEET=$(ls "$SAMPLE_SHEET_DIR"*.csv 2>/dev/null | head -1)
    
    if [ -z "$SAMPLE_SHEET" ]; then
        log_error "No CSV file found in $SAMPLE_SHEET_DIR"
        exit 1
    fi

    check_file_exists "$SAMPLE_SHEET"
    log_info "Found Sample Sheet: $SAMPLE_SHEET"

    # Extract Run ID from filename (format: YYYYMMDD_A*_****_AH****_...)
    RUN_ID=$(basename "$SAMPLE_SHEET" .csv)
    RUN_DIR="${RUN_ID}"
    
    # Extract BCL folder name (first part before .csv)
    BCL_FOLDER="${RUN_ID}"

    log_info "Run ID: $RUN_ID"

    # ========================================================================
    # STEP I: Copy Raw Sample Sheet to BCL folder
    # ========================================================================
    log_info "Step I: Copying Raw Sample Sheet to BCL folder..."

    BCL_TARGET_PATH="${BCL_BASE_PATH}/${BCL_FOLDER}"
    
    if [ ! -d "$BCL_TARGET_PATH" ]; then
        log_info "BCL folder not found at: $BCL_TARGET_PATH"
        log_info "Creating directory structure..."
        mkdir -p "$BCL_TARGET_PATH"
    fi

    cp "$SAMPLE_SHEET" "$BCL_TARGET_PATH/SampleSheet.csv"
    log_info "Sample Sheet copied to: $BCL_TARGET_PATH/SampleSheet.csv"

    # ========================================================================
    # STEP II: Process and modify Sample Sheet
    # ========================================================================
    log_info "Step II: Processing and modifying Sample Sheet..."

    # Create temporary file for processing
    TEMP_SHEET=$(mktemp)
    trap "rm -f $TEMP_SHEET" EXIT

    # Remove [Header], [Reads], and [Settings] sections
    # Keep only [Data] section and remove the [Data] header
    awk '
    BEGIN { in_data = 0 }
    /^\[Data\]/ { in_data = 1; next }
    /^\[Header\]/ || /^\[Reads\]/ || /^\[Settings\]/ { in_data = 0; next }
    in_data && NF > 0 { print }
    ' "$SAMPLE_SHEET" > "$TEMP_SHEET"

    # Process columns
    # 1. Replace Species with Lane (*)
    # 2. Remove Sample_Name
    # 3. Rename Sample_ID → Sample
    # 4. Rename index10X → Index
    # 5. Remove quotes
    
    awk -F',' -v OFS=',' '
    NR == 1 {
        # Process header
        new_header = ""
        for (i = 1; i <= NF; i++) {
            col = $i
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", col)
            
            if (col == "Sample_ID") {
                col = "Sample"
            } else if (col == "index10X") {
                col = "Index"
            } else if (col == "Species") {
                col = "Lane"
            } else if (col == "Sample_Name") {
                continue  # Skip this column
            }
            
            # Remove quotes
            gsub(/"/, "", col)
            
            if (new_header != "") new_header = new_header OFS
            new_header = new_header col
        }
        print new_header
        next
    }
    {
        # Process data rows
        new_row = ""
        for (i = 1; i <= NF; i++) {
            col = $i
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", col)
            
            # Skip Sample_Name column (check header column names)
            if (NR > 1) {
                gsub(/"/, "", col)
                
                # Replace Species values with * for Lane
                if (col == "Species" || (i > 1 && $(1) != "Sample_ID")) {
                    # This is a data cell, replace with *
                    if (col ~ /Species/ || i == NF) {
                        col = "*"
                    }
                }
            }
            
            if (new_row != "") new_row = new_row OFS
            new_row = new_row col
        }
        if (new_row != "") print new_row
    }
    ' "$TEMP_SHEET" > "${TEMP_SHEET}.tmp"
    
    mv "${TEMP_SHEET}.tmp" "$TEMP_SHEET"

    # ========================================================================
    # STEP III: Save modified Sample Sheet
    # ========================================================================
    log_info "Step III: Saving modified Sample Sheet..."

    # Extract base name for output file (YYYYMMDD_A*_****_*)
    OUTPUT_NAME=$(echo "$RUN_ID" | cut -d'_' -f1-4)
    OUTPUT_FILE="Index_mkfastq_${OUTPUT_NAME}.csv"

    # Extract cohort name from run ID (assuming format includes cohort identifier)
    # This may need adjustment based on your actual run ID format
    COHORT=$(echo "$RUN_ID" | grep -oP '_\K[^_]*(?=_)' | head -1)
    COHORT="${COHORT:-default}"

    SHARE_TARGET_PATH="${SHARE_BASE_PATH}/index_mkfastq/${COHORT}"
    mkdir -p "$SHARE_TARGET_PATH"

    cp "$TEMP_SHEET" "$SHARE_TARGET_PATH/$OUTPUT_FILE"
    log_info "Modified Sample Sheet saved to: $SHARE_TARGET_PATH/$OUTPUT_FILE"

    # ========================================================================
    # Completion
    # ========================================================================
    log_info "Pre-processing completed successfully!"
    log_info "Output file: $SHARE_TARGET_PATH/$OUTPUT_FILE"
}

# Execute main function
main "$@"
