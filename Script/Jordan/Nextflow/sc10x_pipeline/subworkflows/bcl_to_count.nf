/*
========================================================================================
    SUB-WORKFLOW : BCL_TO_COUNT
========================================================================================
    Description : Orchestrates BCL -> FASTQ conversion, then alignment / quantification
                  for all samples. Handles per-sample parallelization after mkfastq.
    Steps       :
        1. CELLRANGER_MKFASTQ  (BCL + sample sheet -> FASTQ per sample)
        2. CELLRANGER_COUNT    (FASTQ per sample -> matrices + metrics)
----------------------------------------------------------------------------------------
    Inputs :
        ch_bcl_input : Channel[ tuple(bcl_dir, sample_sheet) ]
    Outputs :
        ch_matrices       : Channel[ tuple(sample_id, filtered_matrix_dir) ]
        ch_metrics        : Channel[ tuple(sample_id, metrics_summary.csv) ]
        ch_web_summaries  : Channel[ tuple(sample_id, web_summary.html) ]
        ch_molecule_info  : Channel[ tuple(sample_id, molecule_info.h5) ]
        ch_versions       : Channel[ versions.yml ]
========================================================================================
*/

include { CELLRANGER_MKFASTQ } from '../modules/cellranger_mkfastq'
include { CELLRANGER_COUNT   } from '../modules/cellranger_count'

workflow BCL_TO_COUNT {

    take:
        ch_bcl_input    // Channel[ tuple(bcl_dir, sample_sheet) ]

    main:

        ch_versions = Channel.empty()

        // -----------------------------------------------------------------------
        // STEP 1: BCL -> FASTQ conversion
        // -----------------------------------------------------------------------
        // Add a run identifier based on BCL directory name
        ch_mkfastq_input = ch_bcl_input.map { bcl_dir, sample_sheet ->
            def run_id = bcl_dir.name.replaceAll(/[^a-zA-Z0-9_-]/, '_')
            return tuple(run_id, bcl_dir, sample_sheet)
        }

        CELLRANGER_MKFASTQ(ch_mkfastq_input)

        ch_versions = ch_versions.mix(CELLRANGER_MKFASTQ.out.versions.first())

        // -----------------------------------------------------------------------
        // STEP 2: Detect FASTQ sample directories after mkfastq
        //
        // mkfastq output structure:
        //   fastq_output/
        //     Sample_A/   <- one directory per sample
        //       Sample_A_S1_L001_R1_001.fastq.gz
        //       Sample_A_S1_L001_R2_001.fastq.gz
        //     Sample_B/
        //       ...
        //
        // Expand the channel to produce one tuple (sample_id, fastq_dir)
        // per detected sample in fastq_output/.
        // -----------------------------------------------------------------------
        ch_count_input = CELLRANGER_MKFASTQ.out.fastqs
            .flatMap { run_id, fastq_dirs ->
                // fastq_dirs can be a Path or a list of Paths
                def dirs = fastq_dirs instanceof List ? fastq_dirs : [fastq_dirs]
                dirs
                    .findAll { it.isDirectory() }
                    .collect { dir ->
                        def sample_id = dir.name
                        log.info "  Sample detected after mkfastq: ${sample_id}"
                        return tuple(sample_id, dir)
                    }
            }

        // Ensure at least one sample was found
        ch_count_input.ifEmpty {
            error "ERROR (BCL_TO_COUNT): No sample directory found after mkfastq. " +
                "Check the sample sheet and BCL data."
        }

        // Shared genome reference for all samples
        ch_genome = Channel.value(file(params.genome_reference, checkIfExists: true))

        // -----------------------------------------------------------------------
        // STEP 3: Alignment and quantification (one parallel job per sample)
        // -----------------------------------------------------------------------
        CELLRANGER_COUNT(ch_count_input, ch_genome)

        ch_versions = ch_versions.mix(CELLRANGER_COUNT.out.versions.first())

    emit:
        matrices      = CELLRANGER_COUNT.out.matrices
        metrics       = CELLRANGER_COUNT.out.metrics
        web_summaries = CELLRANGER_COUNT.out.web_summaries
        molecule_info = CELLRANGER_COUNT.out.molecule_info
        versions      = ch_versions
}
