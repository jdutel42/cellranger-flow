#!/bin/bash
#SBATCH --job-name=nextflow_cellranger
#SBATCH --partition=phoenix
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/home/dutel/Log/nextflow_controller_%j.out
#SBATCH --error=/home/dutel/Log/nextflow_controller_%j.err

cd /labos/UGM/home/dutel/cellranger-flow
module load devel/python/Miniconda3-py39_4.10.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nextflow_env
nextflow run main.nf -c assets/nextflow.config -params-file assets/params.json -profile slurm -resume
