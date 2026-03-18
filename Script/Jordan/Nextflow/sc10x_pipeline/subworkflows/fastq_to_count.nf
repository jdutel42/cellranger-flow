/*
========================================================================================
    SOUS-WORKFLOW : FASTQ_TO_COUNT
========================================================================================
    Description : Orchestre l'alignement / quantification directement à partir de FASTQ
                  pré-existants. Chaque sample est traité en parallèle.
    Étapes      :
        1. CELLRANGER_COUNT (FASTQ par sample → matrices + métriques)
----------------------------------------------------------------------------------------
    Convention attendue pour le dossier d'entrée :
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
        // Validation et logging des samples détectés
        // -----------------------------------------------------------------------
        ch_validated_fastqs = ch_fastq_dirs.map { sample_id, fastq_dir ->

            // Vérification de l'existence du dossier
            if (!fastq_dir.exists()) {
                error "ERREUR (FASTQ_TO_COUNT): Le dossier FASTQ du sample '${sample_id}' n'existe pas : ${fastq_dir}"
            }

            // Vérification qu'au moins un FASTQ R1 et R2 est présent
            def r1_files = fastq_dir.list().findAll { it =~ /_R1_.*\\.fastq\\.gz$/ }
            def r2_files = fastq_dir.list().findAll { it =~ /_R2_.*\\.fastq\\.gz$/ }

            if (r1_files.isEmpty()) {
                log.warn "ATTENTION: Aucun fichier R1 trouvé pour le sample '${sample_id}' dans ${fastq_dir}. Sample ignoré."
                return null
            }
            if (r2_files.isEmpty()) {
                log.warn "ATTENTION: Aucun fichier R2 trouvé pour le sample '${sample_id}' dans ${fastq_dir}. Sample ignoré."
                return null
            }

            log.info "  Sample FASTQ validé : ${sample_id} (${r1_files.size()} fichiers R1, ${r2_files.size()} fichiers R2)"
            return tuple(sample_id, fastq_dir)
        }
        .filter { it != null }

        // Référence génomique partagée entre tous les samples
        ch_genome = Channel.value(file(params.genome_reference, checkIfExists: true))

        // -----------------------------------------------------------------------
        // ÉTAPE 1 : Alignement et quantification (un job par sample en parallèle)
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
