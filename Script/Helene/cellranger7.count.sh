#!/bin/bash
# Cellranger multi script for multiple samples
# One job by sample
# steps : Configuration file, slurm script, job  

########################################################
#           Initialization                             #  
########################################################

# sample are called batch[number]_[CD45plus | Tcell]
# for each batch, two samples : CD45plus and Tcell 

# Input sample name
batch="BM_CREneg BM_CREpos SP_CREpos SP_CREneg" # list of batch number 

# Input
#path_ref=/home/save/hdaunes/cellranger_reference_human/refdata-gex-GRCh38-2020-A #human
path_ref=/home/hdaunes/save/reference_souris/refdata-gex-mm10-2020-A #mouse
path_fastq=/home/hdaunes/work/cellranger_mkfastq_output/HWMK2DSXC

# Output & 
path_libraries=/home/work/hdaunes/libraries # configuration file .csv folder
path_script=/home/hdaunes/work/scr/CR.multi_slurm.jobs # where to save slurm job scripts
output=/home/hdaunes/work/cellranger_count_output # cellranger multi output folder 

#protocol_prefix="RESISTEC" 
protocol_prefix="" 

today=$(date +%Y%m%d) 


########################################################
#           II Writte slurm script                     #
########################################################

# Writte one script by sample eg MIDAS_batch10_CD45plus
for b in ${batch};do

    # Create slurm script 
        script=/home/hdaunes/work/scr/CR.multi_slurm.jobs/CR.count_${b}_${today}.slurm
        test -f $script && rm -r $script # delete slurm script if it exist 

        # Slurm headers
        echo -e "#!/bin/bash \
           \n#SBATCH --partition=workq \
           \n#SBATCH --mail-user=helene.daunes@inserm.fr \
           \n#SBATCH --mail-type=END,FAIL \
           \n#SBATCH --job-name=CRc${b} \
           \n#SBATCH --ntasks=1 \
           \n#SBATCH --cpus-per-task=8 \
           \n#SBATCH --mem=40gb \
           \n#SBATCH --error=/home/hdaunes/work/logs/$today.cr_count.${b}.err \
           \n#SBATCH --output=/home/hdaunes/work/logs/$today.cr_count.${b}.out\n" >> $script


        # Loading modules
        echo -e "module load bioinfo/CellRanger/7.1.0\n" >> $script
 
        # Main function (cellranger count) 
        echo "cellranger count \\
        --id=${b} \\
        --sample=${b} \\
        --fastqs=${path_fastq} \\
        --transcriptome=${path_ref}" >> $script

       

########################################################
#           III Send Job                               #
########################################################
        
    # go to working (=output) directory
        cd ${output}
        
        sbatch ${script}

done





