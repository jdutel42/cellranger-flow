#!/bin/bash

#SBATCH --job-name=cellranger_mkfastq
#SBATCH --partition=phoenix
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb 
#SBATCH --time=24:00:00
#SBATCH --mail-user=jordan.dutel@inserm.fr 
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/dutel/log/%x_%j.out
#SBATCH --error=/home/dutel/log/%x_%j.err   

# Stop the script if a command fails, if an undefined variable is used, or if a command in a pipeline fails
set -euo pipefail 

##########################################
#                  Env                   #
##########################################

# For subshells to have access to conda activate
source /labos/UGM/dev/envs/miniconda/etc/profile.d/conda.sh 

# This is the name of an env with bcl2fastq necessary for cellranger mkfastq to work
conda activate /labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_

##########################################
#                  Input                 #
##########################################

RUN_ID="260304_A01789_0434_BHJKHTDRX7"
BCL_DIR="/sequenceurs/NovaSeq1/${RUN_ID}"
SAMPLE_SHEET="/home/dutel/samplesheet/Index_mkfastq_${RUN_ID}.csv"

##########################################
#                 Output                 #
##########################################

OUT_DIR="/labos/UGM/Recherche/Bi-spe"

##########################################
#                 Verif                 #
##########################################

# Verify that the BCL directory exist
if [ ! -d "${BCL_DIR}" ]; then
    echo "Error: BCL directory ${BCL_DIR} does not exist."
    exit 1
fi

# Verify that the sample sheet exist
if [ ! -f "${SAMPLE_SHEET}" ]; then
    echo "Error: Sample sheet ${SAMPLE_SHEET} does not exist."
    exit 1
fi

# Verify that the output directory exist, if not create it
if [ ! -d "${OUT_DIR}" ]; then
    echo "Output directory ${OUT_DIR} does not exist. Creating it."
    mkdir -p "${OUT_DIR}"
fi

##########################################
#                Logging                 #
##########################################

cd "$OUT_DIR"

echo "========================================"
echo "Job start:"
echo "Date: $(date)"
echo "User: ${USER}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_JOB_NODELIST}"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Working directory: $(pwd)"
echo "----------------------------------------"

echo "Resources:"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 64G"

echo "----------------------------------------"
echo "Inputs:"
echo "RUN_ID: ${RUN_ID}"
echo "BCL_DIR: ${BCL_DIR}"
echo "SAMPLE_SHEET: ${SAMPLE_SHEET}"

echo "----------------------------------------"
echo "Outputs:"
echo "OUT_DIR: ${OUT_DIR}"

echo "----------------------------------------"
echo "Software:"
echo "Cell Ranger version: $(/labos/UGM/dev/cellranger-7.1.0/bin/cellranger --version | head -n 1)"

echo "----------------------------------------"
echo "Command:"

CMD="/labos/UGM/dev/cellranger-7.1.0/bin/cellranger mkfastq \
  --run=${BCL_DIR} \
  --id=${RUN_ID} \
  --csv=${SAMPLE_SHEET} \
  --output-dir=${OUT_DIR} \
  --localcores=${SLURM_CPUS_PER_TASK} \
  --localmem=64"

echo "$CMD"
echo "========================================"

##########################################
#                   Run                  #
##########################################


echo "MKFASTQ START: $(date)"

# Allow the script to continue even if the command fails, so we can capture the exit code and log it
set +e
# Run the command
$CMD
# Capture the exit code of the command
EXIT_CODE=$?
# Re-enable the option to exit on error
set -e

echo "MKFASTQ END: $(date)"
echo "Exit code: $EXIT_CODE"

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: Cellranger failed"
    exit $EXIT_CODE
fi

echo "JOB END: $(date)"