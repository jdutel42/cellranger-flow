#!/usr/bin/env bash
# ==============================================================================
# Script     : Loop_slurm_job_cellranger_multi_9_GEX_VDJ.sh
# Project    : Routine pipeline for single cell RNA-seq data processing (10x Genomics)
# Description: Submit SLURM jobs for cellranger multi (GEX + VDJ) per sample. 
#              Generates a CSV configuration file, a SLURM script, and submits the job. One job per sample.
# Usage      : bash Loop_slurm_job_cellranger_multi_9_GEX_VDJ.sh [--dry-run] [--force]
#              --dry-run  : generates SLURM scripts without submitting jobs (for testing)
#              --force    : overwrites existing SLURM scripts and configuration files without prompting
# Author     : Jordan Dutel
# Email      : jordan.dutel@inserm.fr
# Created    : 2026_03_19
# ==============================================================================

set -euo pipefail  # Exit on error, undefined variable, and fail on pipe errors
IFS=$'\n\t'        # Secure field separator (newline and tab, not space)

# ==============================================================================
#  GLOBAL VARIABLES
# ==============================================================================

readonly SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly TODAY="$(date +%Y%m%d)"
readonly LOG_FILE="/home/dutel/log/${TODAY}.${SCRIPT_NAME%.sh}.log"

# --- Flags CLI ---
DRY_RUN=false
FORCE=false
SLURM_SCRIPT=""

# ==============================================================================
#  CONFIGURATION
# ==============================================================================

# Protocol prefix for sample naming
readonly PROTOCOL_PREFIX="MIDAS_2"

# Batch numbers to process (e.g., BATCH_IDS=(3 4 5 6 7) or BATCH_IDS=(62))
BATCH_IDS=(74 75 76 77 78 79 80)

# --- References ---
readonly PATH_REF_GEX="/labos/UGM/dev/cellranger-pipe/refdata-gex-GRCh38-2020-A"
readonly PATH_REF_VDJ="/labos/UGM/dev/cellranger-pipe/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0"

# --- FASTQ ---
readonly PATH_FASTQ="/labos/UGM/Recherche/midas/fastq2"
readonly FASTQ_FOLDER_GEX="HCHNTDMX2"
readonly FASTQ_FOLDER_VDJ="HCHNTDMX2"

# --- Output paths ---
readonly PATH_SAMPLE_SHEET="/home/dutel/data/samplesheet"   # Sample sheet folder (template and sample-specific configs)
readonly PATH_SAMPLE_SHEET_TEMPLATE="${PATH_SAMPLE_SHEET}/Template" # Template config folder
readonly PATH_SAMPLE_SHEET_SAMPLE="${PATH_SAMPLE_SHEET}/Sample"     # Sample-specific config folder
readonly PATH_SLURM_SCRIPTS="/home/dutel/script/Classique_pipeline/03_Alignment/Cellranger_multi/Cellranger_7/Slurm_job_cellranger_multi_7_GEX_VDJ"  # SLURM scripts folder
readonly PATH_OUTPUT="/labos/UGM/Recherche/midas/output2"   # CellRanger output folder

# --- CellRanger ---
readonly CELLRANGER_BIN="/labos/UGM/dev/cellranger-9.0.1/bin/cellranger"

# --- SLURM ---
readonly SLURM_PARTITION="phoenix"
readonly SLURM_MAIL="jordan.dutel@inserm.fr"
readonly SLURM_CPUS=16
readonly SLURM_MEM="40gb"

# ==============================================================================
#  FUNCTIONS
# ==============================================================================

# Logging
log() {
    local level="$1"; shift
    local msg="$*"
    local ts
    ts="$(date '+%Y-%m-%d %H:%M:%S')"
    echo "[${ts}] [${level}] ${msg}" | tee -a "${LOG_FILE}"
}

log_info()  { log "INFO " "$@"; }
log_warn()  { log "WARN " "$@"; }
log_error() { log "ERROR" "$@" >&2; }

# Exit with error message
die() {
    log_error "$*"
    exit 1
}

# Display help
usage() {
    grep '^# Usage\|^#              ' "${BASH_SOURCE[0]}" | sed 's/^# //'
    exit 0
}

# Check that required tools are available
check_dependencies() {
    local deps=("sbatch" "${CELLRANGER_BIN}")
    for dep in "${deps[@]}"; do
        if ! command -v "${dep}" &>/dev/null && [[ ! -x "${dep}" ]]; then
            die "Missing dependency : ${dep}"
        fi
    done
    log_info "Dependencies checked."
}

# Check that critical directories exist
check_directories() {
    local dirs=("${PATH_REF_GEX}" "${PATH_REF_VDJ}" "${PATH_FASTQ}" \
                "${PATH_SAMPLE_SHEET}" "${PATH_SLURM_SCRIPTS}" "${PATH_OUTPUT}" \
                "${PATH_SAMPLE_SHEET_TEMPLATE}" "${PATH_SAMPLE_SHEET_SAMPLE}" \
                "${PATH_FASTQ}/${FASTQ_FOLDER_GEX}" "${PATH_FASTQ}/${FASTQ_FOLDER_VDJ}")
    for d in "${dirs[@]}"; do
        if [[ ! -d "${d}" ]]; then
            die "Directory not found : ${d}"
        fi
    done
    log_info "Critical directories checked."
}

# ==============================================================================
#  I. GENERATION OF CONFIGURATION FILES
#  - A template config with common settings (references, etc.)
#  - A sample-specific config for each batch with FASTQ paths
# ==============================================================================

# Generate a common template configuration file if it does not exist
generate_config_template() {
    local conf="${PATH_SAMPLE_SHEET_TEMPLATE}/config_template.csv"

    if [[ -f "${conf}" ]]; then
        log_info "Configuration template already exists : ${conf}"
        return 0
    fi

    log_info "Creating configuration template : ${conf}"
    cat > "${conf}" <<EOF
[gene-expression]
ref,${PATH_REF_GEX}
no-bam,FALSE
no-secondary,FALSE

[vdj]
ref,${PATH_REF_VDJ}

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
EOF
}

# Generate a sample-specific configuration file for a given batch
generate_sample_config() {
    local batch="$1"
    local conf_template="${PATH_SAMPLE_SHEET_TEMPLATE}/config_template.csv"
    local conf_sample="${PATH_SAMPLE_SHEET_SAMPLE}/config_sample_${batch}.csv"

    if [[ -f "${conf_sample}" ]]; then
        log_warn "Sample-specific configuration file already exists for ${batch}, overwriting : ${conf_sample}"
        rm "${conf_sample}"
    fi

    log_info "Generating configuration for ${batch} : ${conf_sample}"
    cp "${conf_template}" "${conf_sample}" # Copy template to sample-specific config before appending sample-specific entries

    # Append sample-specific FASTQ entries to the config file
    printf '%s,%s,any,%s,Gene Expression,\n' \
        "${batch}_GEX" \
        "${PATH_FASTQ}/${FASTQ_FOLDER_GEX}" \
        "${batch}_GEX" >> "${conf_sample}"

    printf '%s,%s,any,%s,VDJ,\n' \
        "${batch}_VDJ" \
        "${PATH_FASTQ}/${FASTQ_FOLDER_VDJ}" \
        "${batch}_VDJ" >> "${conf_sample}"
}

# ==============================================================================
#  II. VERIFICATION OF INPUT FASTQ FILES
#  - Check that FASTQ files for GEX and VDJ exist for the given batch
#  - Uses a secure globbing method to avoid issues with special characters
# ==============================================================================

check_fastq() {
    local batch="$1"
    local gex_pattern="${PATH_FASTQ}/${FASTQ_FOLDER_GEX}/${batch}_GEX"
    local vdj_pattern="${PATH_FASTQ}/${FASTQ_FOLDER_VDJ}/${batch}_VDJ"
    local ok=true

    # Use find with -name to check for the presence of FASTQ files matching the expected pattern (glob patterns)
    if ! find "${PATH_FASTQ}/${FASTQ_FOLDER_GEX}" \
            -maxdepth 1 -name "${batch}_GEX*.fastq.gz" | grep -q .; then 
        log_error "FASTQ GEX not found : ${gex_pattern}*.fastq.gz"
        ok=false
    fi

    if ! find "${PATH_FASTQ}/${FASTQ_FOLDER_VDJ}" \
            -maxdepth 1 -name "${batch}_VDJ*.fastq.gz" | grep -q .; then
        log_error "FASTQ VDJ not found : ${vdj_pattern}*.fastq.gz"
        ok=false
    fi

    [[ "${ok}" == true ]]
}

# ==============================================================================
#  III. GENERATION OF SLURM SCRIPTS
#  - For each batch, generate a SLURM script that runs cellranger multi 
#    with the appropriate configuration file
# ==============================================================================

generate_slurm_script() {
    local batch="$1"
    local conf_sample="${PATH_SAMPLE_SHEET_SAMPLE}/config_sample_${batch}.csv"
    # local script="${PATH_SLURM_SCRIPTS}/Slurm_job_cellranger_multi_9_GEX_VDJ_${batch}_${TODAY}.sh"
    SLURM_SCRIPT="${PATH_SLURM_SCRIPTS}/Slurm_job_cellranger_multi_9_GEX_VDJ_${batch}_${TODAY}.sh"


    # Remove existing script if it exists and --force is specified
    if [[ "${FORCE}" == true ]] && [[ -f "${SLURM_SCRIPT}" ]]; then
        log_warn "Script SLURM déjà présent pour ${batch}, suppression forcée : ${SLURM_SCRIPT}"
        rm "${SLURM_SCRIPT}"
    fi

    log_info "SLURM script generation : ${SLURM_SCRIPT}"

    cat > "${SLURM_SCRIPT}" <<EOF
#!/usr/bin/env bash
# ---------------------------------------------------------------
# Job SLURM    : cellranger multi - ${batch}
# Generated by : ${SCRIPT_NAME}
# Date         : ${TODAY}
# ---------------------------------------------------------------
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --mail-user=${SLURM_MAIL}
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cellranger_multi_${batch}
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --error=/home/dutel/log/${TODAY}.Slurm_job_cellranger_multi_9_GEX_VDJ_${batch}.%j.err
#SBATCH --output=/home/dutel/log/${TODAY}.Slurm_job_cellranger_multi_9_GEX_VDJ_${batch}.%j.out

set -euo pipefail

echo "[INFO] Starting cellranger multi slurm job for ${batch}"
echo "[INFO] Node  : \$(hostname)"
echo "[INFO] Date  : \$(date)"

cd "${PATH_OUTPUT}"

${CELLRANGER_BIN} multi \\
    --id="${batch}" \\
    --csv="${conf_sample}" \\
    --localcores=${SLURM_CPUS} \\
    --localmem=38

echo "[INFO] cellranger multi completed for ${batch} — \$(date)"
EOF

    chmod 750 "${SLURM_SCRIPT}"
    echo "${SLURM_SCRIPT}"
}

# ==============================================================================
#  IV. JOB SUBMISSION
#  - Submit the generated SLURM script for each batch
#  - In dry-run mode, only log the intended submission without actually submitting
# ==============================================================================

submit_job() {
    local script="$1"
    local batch="$2"

    if [[ "${DRY_RUN}" == true ]]; then
        log_info "[DRY-RUN] Job not submitted : ${script}"
        return 0
    fi

    local job_id
    job_id="$(sbatch "${script}" | awk '{print $NF}')"
    log_info "Job submitted for ${batch} — SLURM JobID : ${job_id}"
}

# ==============================================================================
#  PARSING CLI ARGUMENTS
#  - Supports --dry-run and --force flags
#  - Displays usage information with -h or --help
# ==============================================================================

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --dry-run) DRY_RUN=false ;;
            --force)   FORCE=false ;;
            -h|--help) usage ;;
            *) die "Unknown argument : $1" ;;
        esac
        shift
    done
}

# ==============================================================================
#  PRIMARY FUNCTION
#  - Orchestrates the entire workflow: configuration generation, verification, script generation, and job submission
#  - Logs the progress and summary of the operations
#  - Exits with an error code if any batch fails to process
# ==============================================================================

main() {
    parse_args "$@"

    # Initialize log file
    mkdir -p "$(dirname "${LOG_FILE}")"

    log_info "========================================================"
    log_info "Démarrage de ${SCRIPT_NAME}"
    log_info "Mode dry-run : ${DRY_RUN} | Force : ${FORCE}"
    log_info "Batches      :$(IFS=','; echo "${BATCH_IDS[*]}")"
    log_info "========================================================"

    check_dependencies
    check_directories

    # Step 1: Generate template configuration if it does not exist
    generate_config_template

    local success=0
    local skipped=0
    local failed=0

    for batch_id in "${BATCH_IDS[@]}"; do

        local batch="${PROTOCOL_PREFIX}_batch${batch_id}"
        log_info "--- Processing ${batch} ---"

        # Step 2: Check output directory
        if [[ -d "${PATH_OUTPUT}/${batch}" ]]; then
            if [[ "${FORCE}" == true ]]; then
                log_warn "Existing output directory found, forcing overwrite : ${PATH_OUTPUT}/${batch}"
                rm -rf "${PATH_OUTPUT}/${batch}"
            else
                log_error "Output directory already exists : ${PATH_OUTPUT}/${batch}"
                log_error "Use --force to overwrite or delete it manually."
                (( failed++ )) || true
                continue
            fi
        fi

        # Step 2: Check FASTQ files
        if ! check_fastq "${batch}"; then
            log_error "FASTQ files missing for ${batch}, batch skipped."
            (( failed++ )) || true
            continue
        fi

        # Step 3: Generate sample-specific configuration
        generate_sample_config "${batch}"

        # Step 4: Generate SLURM script
        generate_slurm_script "${batch}"

        # Step 5: Submit job
        submit_job "${SLURM_SCRIPT}" "${batch}"
        (( success++ )) || true

    done

    # Summary
    log_info "========================================================"
    log_info "Summary : ${success} submitted | ${skipped} skipped | ${failed} failed"
    log_info "Complete log : ${LOG_FILE}"
    log_info "========================================================"

    [[ "${failed}" -eq 0 ]] || exit 1
}

main "$@"