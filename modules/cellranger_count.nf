/*
========================================================================================
    MODULE : CELLRANGER_COUNT
========================================================================================
    Description : Aligns 10x FASTQ files and quantifies per-cell gene expression
                  using cellranger count. Handles samples in parallel.
    Outil       : Cell Ranger 7.2.0
    Container   : nfcore/cellranger:7.2.0
----------------------------------------------------------------------------------------
    Inputs :
        tuple val(sample_id), path(fastq_dir)   → un tuple par sample
        path(genome_reference)                   -> shared across samples
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

    // Publish key outputs to output_dir
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
        // Optional arguments
        def args           = task.ext.args  ?: ''
        def force_cells    = params.force_cells    ? "--force-cells=${params.force_cells}" : ''
        def expect_cells   = params.expect_cells   ? "--expect-cells=${params.expect_cells}" : ''
        def include_introns = params.include_introns ? "--include-introns=true" : "--include-introns=false"

        """
        # -----------------------------------------------------------------------
        # Input validation
        # -----------------------------------------------------------------------
        if [ ! -d "${fastq_dir}" ]; then
            echo "ERROR: FASTQ directory for sample ${sample_id} does not exist: ${fastq_dir}" >&2
            exit 1
        fi

        FASTQ_COUNT=\$(find "${fastq_dir}" -name "*.fastq.gz" | wc -l)
        if [ "\$FASTQ_COUNT" -eq 0 ]; then
            echo "ERROR: No FASTQ files found for sample ${sample_id} in ${fastq_dir}" >&2
            exit 1
        fi

        if [ ! -d "${genome_reference}" ]; then
            echo "ERROR: Genome reference directory does not exist: ${genome_reference}" >&2
            exit 1
        fi

        # -----------------------------------------------------------------------
        # Create log directories
        # -----------------------------------------------------------------------
        mkdir -p logs

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting cellranger count - sample: ${sample_id}" \\
            | tee logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FASTQ dir  : ${fastq_dir}"           | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Genome     : ${genome_reference}"     | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] CPUs       : ${task.cpus}"            | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] RAM        : ${task.memory.toGiga()}G" | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Chemistry  : ${params.chemistry}"     | tee -a logs/count_${sample_id}.log
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Introns    : ${params.include_introns}" | tee -a logs/count_${sample_id}.log

        # -----------------------------------------------------------------------
        # Run cellranger count
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
        # Verify success and required output files
        # -----------------------------------------------------------------------
        if [ \$EXIT_CODE -ne 0 ]; then
            echo "[ERROR] cellranger count failed (code \$EXIT_CODE) for ${sample_id}." | tee -a logs/count_${sample_id}.log
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
                echo "[ERROR] Missing expected output file: \$output_file" | tee -a logs/count_${sample_id}.log
                exit 1
            fi
        done

        # Summarize quantification from the metrics CSV
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] QC metrics for sample ${sample_id}:" | tee -a logs/count_${sample_id}.log
        cat "${sample_id}/outs/metrics_summary.csv" | head -5 | tee -a logs/count_${sample_id}.log

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] cellranger count completed successfully for ${sample_id}." \\
            | tee -a logs/count_${sample_id}.log

        # -----------------------------------------------------------------------
        # Record tool version
        # -----------------------------------------------------------------------
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellranger: \$(cellranger --version 2>&1 | grep -oP 'cellranger-\\K[0-9.]+')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p ${sample_id}/outs/filtered_feature_bc_matrix logs

        # Dummy output files for stub mode
        touch ${sample_id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz
        touch ${sample_id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
        touch ${sample_id}/outs/filtered_feature_bc_matrix/features.tsv.gz
        touch ${sample_id}/outs/molecule_info.h5
        touch ${sample_id}/outs/web_summary.html

        # Dummy metrics_summary.csv
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
