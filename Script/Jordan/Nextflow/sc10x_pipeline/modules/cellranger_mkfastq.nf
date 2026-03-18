/*
========================================================================================
    MODULE : CELLRANGER_MKFASTQ
========================================================================================
    Description : Convertit un dossier BCL en fichiers FASTQ via cellranger mkfastq.
    Outil       : Cell Ranger 7.2.0
    Container   : nfcore/cellranger:7.2.0
----------------------------------------------------------------------------------------
    Inputs :
        tuple val(run_id), path(bcl_dir), path(sample_sheet)
    Outputs :
        tuple val(run_id), path("fastq_output/*")   → ch_fastqs
        path "versions.yml"                          → ch_versions
        path "logs/*.log"                            → ch_logs
========================================================================================
*/

process CELLRANGER_MKFASTQ {

    tag "mkfastq | run: ${run_id}"
    label 'process_high'

    container 'nfcore/cellranger:7.2.0'

    publishDir (
        path    : "${params.output_dir}/mkfastq/${run_id}",
        mode    : 'copy',
        pattern : "fastq_output/**",
        saveAs  : { filename -> filename }
    )
    publishDir (
        path    : "${params.output_dir}/logs/mkfastq",
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    input:
        tuple val(run_id), path(bcl_dir), path(sample_sheet)

    output:
        tuple val(run_id), path("fastq_output/*"), emit: fastqs
        path "versions.yml",                       emit: versions
        path "logs/*.log",                         emit: logs

    when:
        task.ext.when == null || task.ext.when

    script:
        // Construction des arguments optionnels
        def args = task.ext.args ?: ''

        """
        # -----------------------------------------------------------------------
        # Validation des entrées
        # -----------------------------------------------------------------------
        if [ ! -d "${bcl_dir}" ]; then
            echo "ERREUR: Le dossier BCL n'existe pas : ${bcl_dir}" >&2
            exit 1
        fi

        if [ ! -f "${sample_sheet}" ]; then
            echo "ERREUR: La sample sheet n'existe pas : ${sample_sheet}" >&2
            exit 1
        fi

        # -----------------------------------------------------------------------
        # Création des dossiers de sortie et de logs
        # -----------------------------------------------------------------------
        mkdir -p fastq_output logs

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Démarrage cellranger mkfastq — run_id: ${run_id}" \\
            | tee logs/mkfastq_${run_id}.log

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] BCL dir    : ${bcl_dir}"   | tee -a logs/mkfastq_${run_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample sheet: ${sample_sheet}" | tee -a logs/mkfastq_${run_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] CPUs alloués: ${task.cpus}" | tee -a logs/mkfastq_${run_id}.log

        # -----------------------------------------------------------------------
        # Exécution de cellranger mkfastq
        # -----------------------------------------------------------------------
        cellranger mkfastq \\
            --id="${run_id}_mkfastq" \\
            --run="${bcl_dir}" \\
            --csv="${sample_sheet}" \\
            --output-dir=fastq_output \\
            --localcores=${task.cpus} \\
            --localmem=${task.memory.toGiga()} \\
            --delete-undetermined \\
            --ignore-dual-index-flowcells \\
            ${args} \\
            2>&1 | tee -a logs/mkfastq_${run_id}.log

        EXIT_CODE=\${PIPESTATUS[0]}

        # -----------------------------------------------------------------------
        # Vérification du succès
        # -----------------------------------------------------------------------
        if [ \$EXIT_CODE -ne 0 ]; then
            echo "[ERREUR] cellranger mkfastq a échoué avec le code \$EXIT_CODE." | tee -a logs/mkfastq_${run_id}.log
            exit \$EXIT_CODE
        fi

        # Vérification qu'au moins un FASTQ a été généré
        FASTQ_COUNT=\$(find fastq_output -name "*.fastq.gz" | wc -l)
        if [ "\$FASTQ_COUNT" -eq 0 ]; then
            echo "[ERREUR] Aucun fichier FASTQ généré dans fastq_output/" | tee -a logs/mkfastq_${run_id}.log
            exit 1
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] mkfastq terminé — \$FASTQ_COUNT fichiers FASTQ générés." \\
            | tee -a logs/mkfastq_${run_id}.log

        # -----------------------------------------------------------------------
        # Enregistrement de la version de l'outil
        # -----------------------------------------------------------------------
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p fastq_output/Sample_A fastq_output/Sample_B logs

        # Fichiers FASTQ factices pour le mode stub
        touch fastq_output/Sample_A/Sample_A_S1_L001_R1_001.fastq.gz
        touch fastq_output/Sample_A/Sample_A_S1_L001_R2_001.fastq.gz
        touch fastq_output/Sample_B/Sample_B_S2_L001_R1_001.fastq.gz
        touch fastq_output/Sample_B/Sample_B_S2_L001_R2_001.fastq.gz
        touch logs/mkfastq_stub.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: "7.2.0"
        END_VERSIONS
        """
}
