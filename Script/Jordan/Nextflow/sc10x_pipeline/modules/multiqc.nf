/*
========================================================================================
    MODULE : MULTIQC
========================================================================================
    Description : Aggregates all Cell Ranger QC reports (web_summaries, metrics)
                  into a single MultiQC HTML report.
    Outil       : MultiQC 1.21
    Container   : quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0
----------------------------------------------------------------------------------------
    Inputs :
        path(qc_files)       -> List of all QC files to aggregate (collected)
    Outputs :
        path "multiqc_report.html"   → ch_report
        path "multiqc_data/"         → ch_data
        path "versions.yml"          → ch_versions
========================================================================================
*/

process MULTIQC {

    tag "multiqc | global QC aggregation"
    label 'process_medium'

    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    publishDir (
        path    : "${params.output_dir}/multiqc",
        mode    : 'copy',
        saveAs  : { filename -> filename }
    )

    input:
        path qc_files   // Collection of all QC files (metrics_summary.csv, web_summary.html, logs)

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
        # Input checks
        # -----------------------------------------------------------------------
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting MultiQC"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Received QC files:"
        ls -la .

        FILE_COUNT=\$(ls -1 | wc -l)
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] \$FILE_COUNT file(s) to analyze."

        if [ "\$FILE_COUNT" -eq 0 ]; then
            echo "[ERROR] No QC files provided to MultiQC." >&2
            exit 1
        fi

        # -----------------------------------------------------------------------
        # Run MultiQC
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
        # Verify successful execution
        # -----------------------------------------------------------------------
        if [ \$EXIT_CODE -ne 0 ]; then
            echo "[ERROR] MultiQC failed with exit code \$EXIT_CODE." >&2
            exit \$EXIT_CODE
        fi

        if [ ! -f "multiqc_report.html" ]; then
            echo "[ERROR] MultiQC HTML report was not generated." >&2
            exit 1
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] MultiQC completed. Report: multiqc_report.html"

        # -----------------------------------------------------------------------
        # Record tool version
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
