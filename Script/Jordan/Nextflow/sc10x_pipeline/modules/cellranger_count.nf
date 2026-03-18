/*
========================================================================================
    MODULE : CELLRANGER_COUNT
========================================================================================
    Description : Aligne les FASTQ 10x et quantifie l'expression génique par cellule
                  via cellranger count. Gère les samples en parallèle.
    Outil       : Cell Ranger 7.2.0
    Container   : nfcore/cellranger:7.2.0
----------------------------------------------------------------------------------------
    Inputs :
        tuple val(sample_id), path(fastq_dir)   → un tuple par sample
        path(genome_reference)                   → partagé entre samples
    Outputs :
        tuple val(sample_id), path("${sample_id}/outs/filtered_feature_bc_matrix/")
                                                 → ch_matrices
        tuple val(sample_id), path("${sample_id}/outs/metrics_summary.csv")
                                                 → ch_metrics
        tuple val(sample_id), path("${sample_id}/outs/web_summary.html")
                                                 → ch_web_summaries
        tuple val(sample_id), path("${sample_id}/outs/molecule_info.h5")
                                                 → ch_molecule_info (pour cellranger aggr)
        path "versions.yml"                      → ch_versions
        path "logs/*.log"                        → ch_logs
========================================================================================
*/

process CELLRANGER_COUNT {

    tag "count | sample: ${sample_id}"
    label 'process_high'

    container 'nfcore/cellranger:7.2.0'

    // Publier les sorties essentielles dans output_dir
    publishDir (
        path    : "${params.output_dir}/cellranger_count",
        mode    : 'copy',
        pattern : "${sample_id}/outs/**",
        saveAs  : { filename -> filename }
    )
    publishDir (
        path    : "${params.output_dir}/logs/cellranger_count",
        mode    : 'copy',
        pattern : "logs/*.log"
    )

    input:
        tuple val(sample_id), path(fastq_dir)
        path genome_reference

    output:
        tuple val(sample_id), path("${sample_id}/outs/filtered_feature_bc_matrix/"), emit: matrices
        tuple val(sample_id), path("${sample_id}/outs/metrics_summary.csv"),         emit: metrics
        tuple val(sample_id), path("${sample_id}/outs/web_summary.html"),            emit: web_summaries
        tuple val(sample_id), path("${sample_id}/outs/molecule_info.h5"),            emit: molecule_info
        path "versions.yml",                                                          emit: versions
        path "logs/*.log",                                                            emit: logs

    when:
        task.ext.when == null || task.ext.when

    script:
        // Options optionnelles
        def args           = task.ext.args  ?: ''
        def force_cells    = params.force_cells    ? "--force-cells=${params.force_cells}" : ''
        def expect_cells   = params.expect_cells   ? "--expect-cells=${params.expect_cells}" : ''
        def include_introns = params.include_introns ? "--include-introns=true" : "--include-introns=false"

        """
        # -----------------------------------------------------------------------
        # Validation des entrées
        # -----------------------------------------------------------------------
        if [ ! -d "${fastq_dir}" ]; then
            echo "ERREUR: Le dossier FASTQ du sample ${sample_id} n'existe pas : ${fastq_dir}" >&2
            exit 1
        fi

        FASTQ_COUNT=\$(find "${fastq_dir}" -name "*.fastq.gz" | wc -l)
        if [ "\$FASTQ_COUNT" -eq 0 ]; then
            echo "ERREUR: Aucun fichier FASTQ trouvé pour le sample ${sample_id} dans ${fastq_dir}" >&2
            exit 1
        fi

        if [ ! -d "${genome_reference}" ]; then
            echo "ERREUR: Le dossier de référence génomique n'existe pas : ${genome_reference}" >&2
            exit 1
        fi

        # -----------------------------------------------------------------------
        # Création des dossiers de logs
        # -----------------------------------------------------------------------
        mkdir -p logs

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Démarrage cellranger count — sample: ${sample_id}" \\
            | tee logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FASTQ dir  : ${fastq_dir}"           | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Génome     : ${genome_reference}"     | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] CPUs       : ${task.cpus}"            | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] RAM        : ${task.memory.toGiga()}G" | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Chemistry  : ${params.chemistry}"     | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Introns    : ${params.include_introns}" | tee -a logs/count_${sample_id}.log

        # -----------------------------------------------------------------------
        # Exécution de cellranger count
        # -----------------------------------------------------------------------
        cellranger count \\
            --id="${sample_id}" \\
            --fastqs="${fastq_dir}" \\
            --sample="${sample_id}" \\
            --transcriptome="${genome_reference}" \\
            --chemistry="${params.chemistry}" \\
            --localcores=${task.cpus} \\
            --localmem=${task.memory.toGiga()} \\
            ${expect_cells} \\
            ${force_cells} \\
            ${include_introns} \\
            ${args} \\
            2>&1 | tee -a logs/count_${sample_id}.log

        EXIT_CODE=\${PIPESTATUS[0]}

        # -----------------------------------------------------------------------
        # Vérification du succès et des fichiers de sortie attendus
        # -----------------------------------------------------------------------
        if [ \$EXIT_CODE -ne 0 ]; then
            echo "[ERREUR] cellranger count a échoué (code \$EXIT_CODE) pour ${sample_id}." | tee -a logs/count_${sample_id}.log
            exit \$EXIT_CODE
        fi

        REQUIRED_OUTPUTS=(
            "${sample_id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
            "${sample_id}/outs/metrics_summary.csv"
            "${sample_id}/outs/web_summary.html"
            "${sample_id}/outs/molecule_info.h5"
        )

        for output_file in "\${REQUIRED_OUTPUTS[@]}"; do
            if [ ! -f "\$output_file" ] && [ ! -d "\$output_file" ]; then
                echo "[ERREUR] Fichier de sortie attendu manquant : \$output_file" | tee -a logs/count_${sample_id}.log
                exit 1
            fi
        done

        # Résumé de la quantification depuis le CSV des métriques
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Métriques QC du sample ${sample_id} :" | tee -a logs/count_${sample_id}.log
        cat "${sample_id}/outs/metrics_summary.csv" | head -5 | tee -a logs/count_${sample_id}.log

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] cellranger count terminé avec succès pour ${sample_id}." \\
            | tee -a logs/count_${sample_id}.log

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
        mkdir -p ${sample_id}/outs/filtered_feature_bc_matrix logs

        # Fichiers de sortie factices pour le mode stub
        touch ${sample_id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz
        touch ${sample_id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
        touch ${sample_id}/outs/filtered_feature_bc_matrix/features.tsv.gz
        touch ${sample_id}/outs/molecule_info.h5
        touch ${sample_id}/outs/web_summary.html

        # metrics_summary.csv factice
        echo "Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Number of Reads,Valid Barcodes,Sequencing Saturation,Q30 Bases in Barcode,Q30 Bases in RNA Read,Q30 Bases in UMI,Reads Mapped to Genome,Reads Mapped Confidently to Genome,Reads Mapped Confidently to Intergenic Regions,Reads Mapped Confidently to Intronic Regions,Reads Mapped Confidently to Exonic Regions,Reads Mapped Confidently to Transcriptome,Reads Mapped Antisense to Gene,Fraction Reads in Cells,Total Genes Detected,Median UMI Counts per Cell" \\
            > ${sample_id}/outs/metrics_summary.csv
        echo "5000,50000,2500,250000000,98.5%,65.3%,97.1%,96.2%,97.8%,95.4%,91.2%,2.1%,15.3%,73.8%,88.9%,0.5%,85.2%,22000,12000" \\
            >> ${sample_id}/outs/metrics_summary.csv

        touch logs/count_${sample_id}_stub.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: "7.2.0"
        END_VERSIONS
        """
}
