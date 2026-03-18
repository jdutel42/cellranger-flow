#!/bin/bash
# Cellranger multi script for multiple samples
# One job by sample
# steps : Configuration file, slurm script,job  

########################################################
#           Initialization                             #  
########################################################

# sample are called batch[number]_[CD45plus | Tcell]
# for each batch, two samples : CD45plus and Tcell 

# Input sample name
protocol_prefix="MIDAS_2" 
#batch = $(seq -f "${protocol_prefix}_%g" 49 51 )
#batch=$(for i in $(seq 39 47); do echo -n "${protocol_prefix}_batch${i}X "; done) # output_4
batch=$(for i in 39 45 46 47; do echo -n "${protocol_prefix}_batch${i}X "; done) # output_4
batch="MIDAS_2_batch45X"
# Reference
path_ref_gex=/labos/UGM/dev/cellranger-pipe/refdata-gex-GRCh38-2020-A
path_ref_vdj=/labos/UGM/dev/cellranger-pipe/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0

# Input fastq
path_fastq=/labos/UGM/Recherche/midas/fastqX
fastq_folder_gex=H7MNGDMX2
fastq_folder_vdj=${fastq_folder_gex}

# Output & intermediary files
path_libraries=/home/daunes/samplesheet # configuration file .csv folder
path_script=/home/daunes/script/cr_multi_jobs # where to save slurm job scripts
output=/labos/UGM/Recherche/midas/outputX3 # cellranger output folder 

today=$(date +%Y%m%d) 

########################################################
#           I Configuration file                       #
########################################################

# Configuration file template (to be ran 1 time)
conf=${path_libraries}/config_v9.0.1.csv
if ! [ -f ${conf} ]; then
    echo -e "[gene-expression] 
    ref,${path_ref_gex}
    create-bam,TRUE
    no-secondary,FALSE \n\n ">>$conf

    echo -e "[vdj] 
    ref,${path_ref_vdj} \n\n">> $conf


    echo "[libraries]">>$conf
    echo "fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate" >>$conf
fi



for b in ${batch};do
    # copy and rename template according batch number (b) and cell type (t). 
    conf2=${path_libraries}/config.${b}.csv
    #echo $conf2
    if [ -f ${conf2} ]; then
         rm -i $conf2 # ask to rm configuration if the file exists
    fi
    # Specify all input library data (sample specific) if configuration file does not exist
    if ! [ -f ${conf2} ]; then
        cp ${conf} ${conf2} 
        echo "${b}_GEX, ${path_fastq}/${fastq_folder_gex}/$b,any,${b}_GEX, Gene Expression,  ">> $conf2
        echo "${b}_VDJ, ${path_fastq}/${fastq_folder_vdj}/$b,any,${b}_VDJ, VDJ,  ">> $conf2
    fi
    #nano $conf2
done

########################################################
#           II Writte slurm script                     #
########################################################

# Writte one script by sample eg MIDAS_batch10_CD45plus
for b in ${batch};do
    echo $b
    # check for GEX fastq
    #fastqgex=("${path_fastq}/${fastq_folder_gex}/${b}_GEX"*.fastq.gz)
    #if [ ${#fastqgex[@]} -eq 0 ]; then
    #    echo "missing fastq GEX in ${fastq_folder_gex}"
    #    continue
    #fi
    # check for vdj fastq 
    #fastqvdj=("${path_fastq}/${fastq_folder_vdj}/${b}_VDJ"*.fastq.gz)
    #if [ ${#fastqvdj[@]} -eq 0 ]; then
    #    echo "missing fastq VDJ in ${fastq_folder_vdj}"
    #    continue
    #fi
    # check for output folder
    #if [ -d ${output}/${b} ]; then
    #    echo Error output folder already exist
    #    exit -1
    #fi 
    # Create slurm script 
    script=/home/daunes/script/cr_multi_jobs/cr.multi_${b}_${today}.sh
    test -f $script && rm -r $script # delete slurm script if it exist 

    # Slurm headers
    echo -e "#!/bin/bash \
    \n#SBATCH --partition=phoenix \
    \n#SBATCH --mail-user=helene.daunes@inserm.fr \
    \n#SBATCH --mail-type=END,FAIL \
    \n#SBATCH --job-name=CRm${b} \
    \n#SBATCH --cpus-per-task=12 \
    \n#SBATCH --mem=40gb \
    \n#SBATCH --error=/home/daunes/logs/$today.cr_multi.${b}.%j.err \
    \n#SBATCH --output=/home/daunes/logs/$today.cr_multi.${b}.%j.out\n" >> $script

    # Main function (cellranger multi) 
    echo "/labos/UGM/dev/cellranger-9.0.1/bin/cellranger multi \\
    --id=${b} \\
    --csv=${conf2}" >> $script     

########################################################
#           III Send Job                               #
########################################################
   # nano $script 
    # go to working (=output) directory
    cd ${output}
        
    sbatch ${script}
  
done

