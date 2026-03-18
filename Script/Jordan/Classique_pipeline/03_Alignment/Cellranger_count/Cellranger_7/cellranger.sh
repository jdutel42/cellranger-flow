
#!/bin/bash
source /labos/UGM/dev/envs/miniconda/etc/profile.d/conda.sh                 # pour que les subshell aient accès au conda activate           
#conda activate /labos/UGM/dev/envs/shared/178265b579c72c6695d48557d4eadac6_ # c'est le nom d'un env avec bcl2fastq nécessaire à activer pour que cellranger mkfastq fonctionne

#/labos/UGM/dev/cellranger-7.1.0/bin/cellranger multi
batch=$(for i in $(seq 102 102); do echo -n "RESISTEC_batch$i "; done)


# Input
#path_ref=/home/hdaunes/save/reference_souris/refdata-gex-mm10-2020-A
path_ref_gex=/labos/UGM/dev/cellranger-pipe/refdata-gex-GRCh38-2020-A
#path_ref_vdj=/home/hdaunes/save/reference_souris/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0
path_ref_vdj=/labos/UGM/dev/cellranger-pipe/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.0.0
path_fastq=/labos/UGM/Recherche/Resistec/fastq
fastq_folder_gex=H2MG5DMX2
fastq_folder_vdj=H2LMLDMX2 

# Output & 
path_libraries=/home/daunes/samplesheet # configuration file .csv folder
path_script=/home/daunes/script/cr_multi_jobs # where to save slurm job scripts
output=/labos/UGM/Recherche/Resistec/output # cellranger multi output folder 

today=$(date +%Y%m%d) 
for b in $batch;do
   # Error if output folder exist => cellranger will fail 
    if [ -d ${output}/${b} ]; then
        echo Error output folder already exist
        exit -1
    fi 
    # skip batch if vdj fastq not found
    if  [ $(ls ${path_fastq}/${fastq_folder_vdj}/${b}_VDJ_* | wc -l) -eq "0" ];then
        echo Warning $b - VDJ fastq not found in ${fastq_folder_vdj}
        break
    fi
    # skip batch if gex fastq not found
    if ! [ -f ${path_fastq}/${fastq_folder_gex}/${b}_GEX_*.fastq.gz ];then
        echo Warning $b - GEX fastq not found in ${fastq_folder_gex}
        break
    fi 
done