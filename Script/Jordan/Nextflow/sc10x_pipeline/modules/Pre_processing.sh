#!/bin/bash

###############################################################################
# Pre-processing Script for Single Cell Sample Sheet
# Purpose: Prepare and modify SampleSheet for 10X analysis
###############################################################################

set -euo pipefail # Exit on error (set -e), undefined variable (set -u), or error in pipeline (set -o pipefail)

# ============================================================================
# Configuration
# ============================================================================

RAW_SAMPLE_SHEET_FILE_NAME="20260323_scMIDAS2_74_80.csv"
RAW_SAMPLE_SHEET_DIR="/home/dutel/Data/Sample_sheet/Raw"
RAW_SAMPLE_SHEET_FILE_PATH="${RAW_SAMPLE_SHEET_DIR}/${RAW_SAMPLE_SHEET_FILE_NAME}"

RUN_ID="test"
MODIFIED_SAMPLE_SHEET_FILE_NAME="Index_mkfastq_${RUN_ID}.csv"
MODIFIED_SAMPLE_SHEET_DIR="/home/dutel/Data/Sample_sheet/Modified"
MODIFIED_SAMPLE_SHEET_FILE_PATH="${MODIFIED_SAMPLE_SHEET_DIR}/${MODIFIED_SAMPLE_SHEET_FILE_NAME}"

#BCL_BASE_PATH="/media/CRCT13/60To/SingleCell/BclNovaseq"
#SHARE_BASE_PATH="/media/CRCT13/20To"

# ============================================================================
# Functions
# ============================================================================

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

# ============================================================================
# Main Script
# ============================================================================

main() {
    log_info "Starting Sample Sheet preprocessing..."

    # ========================================================================
    # STEP I: Check path and file sample sheet exist
    # ========================================================================
    log_info "Step I: Check path and file sample sheet exist..."

    # Check if Raw Sample Sheet directory exists
    if [ ! -d "$RAW_SAMPLE_SHEET_DIR" ]; then
        log_error "Raw Sample Sheet directory not found: $RAW_SAMPLE_SHEET_DIR"
        exit 1
    fi

    # Check if Raw Sample Sheet file exists
    if [ ! -f "$RAW_SAMPLE_SHEET_FILE_PATH" ]; then
        log_error "$RAW_SAMPLE_SHEET_FILE_NAME file not found in $RAW_SAMPLE_SHEET_DIR"
        exit 1
    fi

    log_info "Found Raw Sample Sheet: $RAW_SAMPLE_SHEET_FILE_PATH"

    # Check if Modified Sample Sheet directory exists
    if [ ! -d "$MODIFIED_SAMPLE_SHEET_DIR" ]; then
        log_error "Modified Sample Sheet directory not found: $MODIFIED_SAMPLE_SHEET_DIR"
        exit 1
    fi

    # Extract Run ID from filename (format: YYYYMMDD_A*_****_AH****_...)
    #RUN_ID=$(basename "$SAMPLE_SHEET" .csv)
    #RUN_DIR="${RUN_ID}"
    
    # Extract BCL folder name (first part before .csv)
    #BCL_FOLDER="${RUN_ID}"

    log_info "Run ID: $RUN_ID"

    # ========================================================================
    # STEP I: Copy Raw Sample Sheet to BCL folder
    # ========================================================================
    #log_info "Step I: Copying Raw Sample Sheet to BCL folder..."

    #BCL_TARGET_PATH="${BCL_BASE_PATH}/${BCL_FOLDER}"
    
    #if [ ! -d "$BCL_TARGET_PATH" ]; then
    #    log_info "BCL folder not found at: $BCL_TARGET_PATH"
    #    log_info "Creating directory structure..."
    #    mkdir -p "$BCL_TARGET_PATH"
    #fi

    #cp "$SAMPLE_SHEET" "$BCL_TARGET_PATH/SampleSheet.csv"
    #log_info "Sample Sheet copied to: $BCL_TARGET_PATH/SampleSheet.csv"

    # ========================================================================
    # STEP II: Process and modify Sample Sheet
    # ========================================================================
    log_info "Step II: Processing and modifying Raw Sample Sheet..."

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
    ' "$RAW_SAMPLE_SHEET_FILE_PATH" > "$TEMP_SHEET"

    # Process columns
    # 1. Replace Species with Lane (*)
    # 2. Remove Sample_Name
    # 3. Rename Sample_ID → Sample
    # 4. Rename index10X → Index
    # 5. Remove quotes
    
    awk -F',' -v OFS=',' '
    NR == 1 {
        # Détection des colonnes
        for (i = 1; i <= NF; i++) {
            col = $i
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", col)
            gsub(/"/, "", col)

            if (col == "Sample_ID") {
                sample_idx = i
                header[i] = "Sample"
            }
            else if (col == "index10X") {
                header[i] = "Index"
            }
            else if (col == "Species") {
                lane_idx = i
                header[i] = "Lane"
            }
            else if (col == "Sample_Name") {
                skip[i] = 1
            }
            else {
                header[i] = col
            }
        }

        # Print header sans colonnes supprimées
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

            col = $i
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", col)
            gsub(/"/, "", col)

            if (i == lane_idx) {
                col = "*"
            }

            if (out != "") out = out OFS
            out = out col
        }
        print out
    }
    ' "$TEMP_SHEET" > "${TEMP_SHEET}.tmp"
    
    mv "${TEMP_SHEET}.tmp" "$TEMP_SHEET"

    # ========================================================================
    # STEP III: Save modified Sample Sheet
    # ========================================================================
    log_info "Step III: Saving modified Sample Sheet..."

    # Extract base name for output file (YYYYMMDD_A*_****_*)
    #OUTPUT_NAME=$(echo "$RUN_ID" | cut -d'_' -f1-4)
    #OUTPUT_FILE="Index_mkfastq_${OUTPUT_NAME}.csv"

    # Extract cohort name from run ID (assuming format includes cohort identifier)
    # This may need adjustment based on your actual run ID format
    #COHORT=$(echo "$RUN_ID" | grep -oP '_\K[^_]*(?=_)' | head -1)
    #COHORT="${COHORT:-default}"

    #SHARE_TARGET_PATH="${SHARE_BASE_PATH}/index_mkfastq/${COHORT}"
    #mkdir -p "$SHARE_TARGET_PATH"

    cp "$TEMP_SHEET" "$MODIFIED_SAMPLE_SHEET_FILE_PATH"
    log_info "Modified Sample Sheet saved to: $MODIFIED_SAMPLE_SHEET_FILE_PATH"

    # ========================================================================
    # Completion
    # ========================================================================
    log_info "Pre-processing completed successfully!"
    log_info "Output file: $MODIFIED_SAMPLE_SHEET_FILE_PATH"
}

# Execute main function
main "$@"
