/*
========================================================================================
    MODULE : MULTIQC
========================================================================================
    Description : Agrège tous les rapports QC de Cell Ranger (web_summaries, metrics)
                  en un rapport MultiQC HTML unique.
    Outil       : MultiQC 1.21
    Container   : quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0
----------------------------------------------------------------------------------------
    Inputs :
        path(qc_files)       → Liste de tous les fichiers QC à agréger (collected)
    Outputs :
        path "multiqc_report.html"   → ch_report
        path "multiqc_data/"         → ch_data
        path "versions.yml"          → ch_versions
========================================================================================
*/

process MULTIQC {

    tag "multiqc | agrégation QC globale"
    label 'process_medium'

    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    publishDir (
        path    : "${params.output_dir}/multiqc",
        mode    : 'copy',
        saveAs  : { filename -> filename }
    )

    input:
        path qc_files   // Collection de tous les fichiers QC (metrics_summary.csv, web_summary.html, logs)

    output:
        path "multiqc_report.html", emit: report
        path "multiqc_data/",       emit: data
        path "versions.yml",        emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args           = task.ext.args ?: ''
        def multiqc_title  = params.multiqc_title  ? "--title \"${params.multiqc_title}\"" : ''
        def multiqc_config = params.multiqc_config && file(params.multiqc_config).exists()
                           ? "--config ${params.multiqc_config}" : ''

        """
        # -----------------------------------------------------------------------
        # Vérification des inputs
        # -----------------------------------------------------------------------
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Démarrage MultiQC"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fichiers QC reçus :"
        ls -la .

        FILE_COUNT=\$(ls -1 | wc -l)
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] \$FILE_COUNT fichier(s) à analyser."

        if [ "\$FILE_COUNT" -eq 0 ]; then
            echo "[ERREUR] Aucun fichier QC fourni à MultiQC." >&2
            exit 1
        fi

        # -----------------------------------------------------------------------
        # Exécution de MultiQC
        # -----------------------------------------------------------------------
        multiqc \\
            ${multiqc_title} \\
            ${multiqc_config} \\
            --force \\
            --verbose \\
            --dirs \\
            --fullnames \\
            --outdir . \\
            ${args} \\
            .

        EXIT_CODE=\$?

        # -----------------------------------------------------------------------
        # Vérification du succès
        # -----------------------------------------------------------------------
        if [ \$EXIT_CODE -ne 0 ]; then
            echo "[ERREUR] MultiQC a échoué avec le code \$EXIT_CODE." >&2
            exit \$EXIT_CODE
        fi

        if [ ! -f "multiqc_report.html" ]; then
            echo "[ERREUR] Le rapport MultiQC HTML n'a pas été généré." >&2
            exit 1
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] MultiQC terminé. Rapport : multiqc_report.html"

        # -----------------------------------------------------------------------
        # Enregistrement de la version de l'outil
        # -----------------------------------------------------------------------
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$(multiqc --version 2>&1 | grep -oP 'version \\K[0-9.]+')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p multiqc_data

        cat <<-EOF > multiqc_report.html
        <html><body><h1>MultiQC Report (stub)</h1></body></html>
        EOF

        touch multiqc_data/multiqc_general_stats.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: "1.21"
        END_VERSIONS
        """
}
