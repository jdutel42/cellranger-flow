/*
========================================================================================
    SUB-WORKFLOW : FASTQ_TO_COUNT
========================================================================================
    Description : Orchestrates alignment / quantification directly from existing
                  FASTQ files. Each sample is processed in parallel.
    Steps       :
        1. CELLRANGER_COUNT (FASTQ per sample -> matrices + metrics)
----------------------------------------------------------------------------------------
    Expected input directory convention:
        <input_dir>/
          <sample_id_1>/
            <sample_id_1>_S1_L001_R1_001.fastq.gz
            <sample_id_1>_S1_L001_R2_001.fastq.gz
            ...
          <sample_id_2>/
            ...
    Inputs :
        ch_fastq_dirs : Channel[ tuple(sample_id, fastq_dir) ]
    Outputs :
        ch_matrices       : Channel[ tuple(sample_id, filtered_matrix_dir) ]
        ch_metrics        : Channel[ tuple(sample_id, metrics_summary.csv) ]
        ch_web_summaries  : Channel[ tuple(sample_id, web_summary.html) ]
        ch_molecule_info  : Channel[ tuple(sample_id, molecule_info.h5) ]
        ch_versions       : Channel[ versions.yml ]
========================================================================================
*/

include { CELLRANGER_COUNT } from '../modules/cellranger_count'

workflow FASTQ_TO_COUNT {

    take:
        ch_fastq_dirs   // Channel[ tuple(sample_id, fastq_dir) ]

    main:

        ch_versions = Channel.empty()

        // -----------------------------------------------------------------------
        // Validate and log detected samples
        // -----------------------------------------------------------------------
        ch_validated_fastqs = ch_fastq_dirs.map { sample_id, fastq_dir ->

            // Check directory exists
            if (!fastq_dir.exists()) {
                error "ERROR (FASTQ_TO_COUNT): FASTQ directory for sample '${sample_id}' does not exist: ${fastq_dir}"
            }

            // Ensure at least one R1 and one R2 FASTQ are present
            def r1_files = fastq_dir.list().findAll { it =~ /_R1_.*\\.fastq\\.gz$/ }
            def r2_files = fastq_dir.list().findAll { it =~ /_R2_.*\\.fastq\\.gz$/ }

            if (r1_files.isEmpty()) {
                log.warn "WARNING: No R1 file found for sample '${sample_id}' in ${fastq_dir}. Sample skipped."
                return null
            }
            if (r2_files.isEmpty()) {
                log.warn "WARNING: No R2 file found for sample '${sample_id}' in ${fastq_dir}. Sample skipped."
                return null
            }

            log.info "  FASTQ sample validated: ${sample_id} (${r1_files.size()} R1 files, ${r2_files.size()} R2 files)"
            return tuple(sample_id, fastq_dir)
        }
        .filter { it != null }

        // Shared genome reference for all samples
        ch_genome = Channel.value(file(params.genome_reference, checkIfExists: true))

        // -----------------------------------------------------------------------
        // STEP 1: Alignment and quantification (one parallel job per sample)
        // -----------------------------------------------------------------------
        CELLRANGER_COUNT(ch_validated_fastqs, ch_genome)

        ch_versions = ch_versions.mix(CELLRANGER_COUNT.out.versions.first())

    emit:
        matrices      = CELLRANGER_COUNT.out.matrices
        metrics       = CELLRANGER_COUNT.out.metrics
        web_summaries = CELLRANGER_COUNT.out.web_summaries
        molecule_info = CELLRANGER_COUNT.out.molecule_info
        versions      = ch_versions
}
