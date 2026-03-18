/*
========================================================================================
    SOUS-WORKFLOW : BCL_TO_COUNT
========================================================================================
    Description : Orchestre la conversion BCL → FASTQ puis l'alignement / quantification
                  pour tous les samples. Gère la parallélisation par sample après mkfastq.
    Étapes      :
        1. CELLRANGER_MKFASTQ  (BCL + sample sheet → FASTQ par sample)
        2. CELLRANGER_COUNT    (FASTQ par sample → matrices + métriques)
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
        // ÉTAPE 1 : Conversion BCL → FASTQ
        // -----------------------------------------------------------------------
        // Ajout d'un identifiant de run basé sur le nom du dossier BCL
        ch_mkfastq_input = ch_bcl_input.map { bcl_dir, sample_sheet ->
            def run_id = bcl_dir.name.replaceAll(/[^a-zA-Z0-9_-]/, '_')
            return tuple(run_id, bcl_dir, sample_sheet)
        }

        CELLRANGER_MKFASTQ(ch_mkfastq_input)

        ch_versions = ch_versions.mix(CELLRANGER_MKFASTQ.out.versions.first())

        // -----------------------------------------------------------------------
        // ÉTAPE 2 : Détection des dossiers FASTQ par sample après mkfastq
        //
        // Structure de sortie de mkfastq :
        //   fastq_output/
        //     Sample_A/   ← un dossier par sample
        //       Sample_A_S1_L001_R1_001.fastq.gz
        //       Sample_A_S1_L001_R2_001.fastq.gz
        //     Sample_B/
        //       ...
        //
        // On "explose" le channel pour obtenir un tuple (sample_id, fastq_dir)
        // par sample détecté dans fastq_output/.
        // -----------------------------------------------------------------------
        ch_count_input = CELLRANGER_MKFASTQ.out.fastqs
            .flatMap { run_id, fastq_dirs ->
                // fastq_dirs peut être un Path ou une liste de Paths
                def dirs = fastq_dirs instanceof List ? fastq_dirs : [fastq_dirs]
                dirs
                    .findAll { it.isDirectory() }
                    .collect { dir ->
                        def sample_id = dir.name
                        log.info "  Sample détecté après mkfastq : ${sample_id}"
                        return tuple(sample_id, dir)
                    }
            }

        // Vérification qu'au moins un sample a été trouvé
        ch_count_input.ifEmpty {
            error "ERREUR (BCL_TO_COUNT): Aucun dossier de sample trouvé après mkfastq. " +
                  "Vérifiez la sample sheet et les données BCL."
        }

        // Référence génomique partagée entre tous les samples
        ch_genome = Channel.value(file(params.genome_reference, checkIfExists: true))

        // -----------------------------------------------------------------------
        // ÉTAPE 3 : Alignement et quantification (un job par sample en parallèle)
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
