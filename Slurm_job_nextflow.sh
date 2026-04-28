#!/bin/bash
#SBATCH --job-name=nextflow_cellranger
#SBATCH --partition=phoenix
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=Log/nextflow_controller_%j.out
#SBATCH --error=Log/nextflow_controller_%j.err

cd /labos/UGM/home/dutel/CellRanger_Flow
source /chemin/vers/votre/env_nextflow/bin/activate
nextflow run main.nf -c assets/nextflow.config -params-file assets/params.json -profile slurm -resume