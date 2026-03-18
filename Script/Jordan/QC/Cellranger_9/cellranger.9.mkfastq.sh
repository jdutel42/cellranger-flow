#!/bin/bash 
#SBATCH --partition=phoenix
#SBATCH --mail-user=jordan.dutel@inserm.fr 
#SBATCH --mail-type=END,FAIL 
#SBATCH --job-name=Mkfastq
#SBATCH --mem=64gb 
#SBATCH --error=/home/daunes/logs/cr_mkfastq_%x_%j.err   
#SBATCH --output=/home/daunes/logs/cr_mkfastq_%x_%j.out


########################################################
#           Initialization                             #  
########################################################

# load modules
source /labos/UGM/dev/envs/miniconda/etc/profile.d/conda.sh                 # pour que les subshell aient accès au conda activate           
conda activate /labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_ # c'est le nom d'un env avec bcl2fastq nécessaire à activer pour que cellranger mkfastq fonctionne

# Input 
path_to_bcl=/sequenceurs/NovaSeq1
bcl=260304_A01789_0434_BHJKHTDRX7
path_libraries=/home/daunes/samplesheet
today=$(date +%Y%m%d) 

#Output
path_fastq=/labos/UGM/Recherche/Bi-spe

########################################################
#           Main                                      #  
########################################################

cd ${path_fastq}

# Main function (cellranger multi) 
/labos/UGM/dev/cellranger-9.0.1/bin/cellranger mkfastq --id=$bcl \
    --run=${path_to_bcl}/${bcl} \
    --csv=${path_libraries}/Index_mkfastq_${bcl}.csv \
    --output-dir=${path_fastq} \
    -p $SLURM_CPUS_PER_TASK -r $SLURM_CPUS_PER_TASK -w $SLURM_CPUS_PER_TASK


#$SLURM_CPUS_PER_TASK
#sbatch -c 48 script.sh