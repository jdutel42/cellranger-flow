#!/bin/bash
# Cellranger multi script for multiple samples
# One job by sample
# steps : Configuration file, slurm script,job  

########################################################
#           Initialization                             #  
########################################################

# Input sample name
#batch=$(for i in $(seq 3 2 15); do echo -n "BR22batch$i "; done) # output_4
# batch=$(for tissu in LN T; do for condition in KO WT; do for treatment in IG BIT; do echo -n $tissu )
batch=(
LN_KO_IG
LN_KO_BIT
LN_WT_IG
LN_WT_BIT
T_KO_IG
T_KO_BIT
T_WT_IG
T_WT_BIT
)

#batch=BR22batch1
# Input
path_ref_gex=/home/daunes/reference_souris/refdata-gex-GRCm39-2024-A
#path_ref_gex=/labos/UGM/dev/cellranger-pipe/refdata-gex-GRCh38-2020-A
path_ref_vdj=/home/daunes/reference_souris/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0
#path_ref_vdj=/labos/UGM/dev/cellranger-pipe/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0
path_ref_feature=/home/daunes/samplesheet/feature_ref_AC_totalseqC.csv

path_fastq=/labos/UGM/Recherche/Bi-spe
fastq_folder_gex=HJKHTDRX7
fastq_folder_vdj=${fastq_folder_gex}
fastq_folder_AC=${fastq_folder_gex}

# Output & 
path_libraries=/home/daunes/samplesheet # configuration file .csv folder
path_script=/home/daunes/script/cr_multi_jobs # where to save slurm job scripts
output=/labos/UGM/Recherche/Bi-spe/output_align # cellranger multi output folder 

today=$(date +%Y%m%d) 

########################################################
#           I Configuration file                       #
########################################################

# Configuration file template (to be ran 1 time)
conf=${path_libraries}/config_AC_mouse.csv
if ! [ -f ${conf} ]; then
    echo -e "[gene-expression] 
ref,${path_ref_gex}
create-bam,TRUE
no-secondary,FALSE 
check-library-compatibility,FALSE\n\n">>$conf

echo -e "[feature]
ref,${path_ref_feature} \n\n">> $conf

    echo -e "[vdj] 
ref,${path_ref_vdj} \n\n">> $conf


    echo "[libraries]">>$conf
    echo "fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate" >>$conf
fi



for b in ${batch};do
    # copy and rename template according batch number (b) and cell type (t). 
        
    conf2=${path_libraries}/config.${b}.csv
    if [ -f ${conf2} ]; then
        #echo existing config file found
        rm $conf2
    fi
    # Specify all input library data (sample specific) if configuration file does not exist
    cp ${conf} ${conf2}
    echo "${b}_GEX, ${path_fastq}/${fastq_folder_gex},any,${b}_GEX, Gene Expression,  ">> $conf2
    echo "${b}_VDJ, ${path_fastq}/${fastq_folder_vdj},any,${b}_VDJ, VDJ-T,  ">> $conf2
    echo "${b}_prot, ${path_fastq}/${fastq_folder_AC},any,${b}_prot, Antibody Capture,  ">> $conf2
    #vim $conf2

    echo -e "\n[samples] 
    sample_id,hashtag_ids">> $conf2
    echo "M2,TotalSeqC_Hashtag_C0302">> $conf2
    echo "M3,TotalSeqC_Hashtag_C0303">> $conf2
    echo "M4,TotalSeqC_Hashtag_C0304">> $conf2
    echo "M5,TotalSeqC_Hashtag_C0305">> $conf2
    echo "M6,TotalSeqC_Hashtag_C0306">> $conf2

done

echo conf done 

########################################################
#           II Writte slurm script                     #
########################################################

# Writte one script by sample eg MIDAS_batch10_CD45plus
for b in ${batch};do
    echo start chek input for $b
    # Error if output folder exist => cellranger will fail 
    if [ -d ${output}/${b} ]; then
        echo Error output folder already exist
        exit -1
    fi 
    # skip batch if gex fastq not found
    if  [ $(ls ${path_fastq}/${fastq_folder_gex} | grep ${b}_GEX* | wc -l) -eq "0" ];then
        echo Warning $b - GEX fastq not found in ${fastq_folder_gex}
        rm ${conf2}
        continue
    fi 
    # path de configuration file
    conf2=${path_libraries}/config.${b}.csv
    # Create slurm script 
    script=/home/daunes/script/cr_multi_jobs/cr.multi_${b}_${today}.sh
    test -f $script && rm -r $script # delete slurm script if it exist 
    # Slurm headers
    echo -e "#!/bin/bash \
    \n#SBATCH --partition=phoenix \
    \n#SBATCH --mail-user=jordan.dutel@inserm.fr \
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
        
    # go to working (=output) directory
    echo $b $script
    cd ${output}       
    sbatch ${script}

done