#!/usr/bin/env nextflow

/*
========================================================================================
    sc10x_pipeline : Single-Cell 10x Genomics Pipeline (DSL2)
========================================================================================
    Auteur    : Bioinformatics Pipeline Team
    Version   : 1.0.0
    Licence   : MIT
    Usage     : nextflow run main.nf -profile docker --input_type fastq --input_dir /data/fastq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// ========================================================================================
// IMPORTS DES MODULES
// ========================================================================================

include { CELLRANGER_MKFASTQ  } from './modules/cellranger_mkfastq'
include { CELLRANGER_COUNT    } from './modules/cellranger_count'
include { MULTIQC             } from './modules/multiqc'

// ========================================================================================
// IMPORTS DES SOUS-WORKFLOWS
// ========================================================================================

include { BCL_TO_COUNT        } from './subworkflows/bcl_to_count'
include { FASTQ_TO_COUNT      } from './subworkflows/fastq_to_count'

// ========================================================================================
// AFFICHAGE DU BANDEAU DE DÉMARRAGE
// ========================================================================================

def printHeader() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║         sc10x_pipeline — Single-Cell 10x Genomics DSL2          ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  input_type    : ${params.input_type}
    ║  input_dir     : ${params.input_dir}
    ║  output_dir    : ${params.output_dir}
    ║  genome        : ${params.genome_reference}
    ║  sample_sheet  : ${params.sample_sheet}
    ║  localcores    : ${params.localcores}
    ║  localmemory   : ${params.localmemory} GB
    ╚══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

// ========================================================================================
// VALIDATION DES PARAMÈTRES OBLIGATOIRES
// ========================================================================================

def validateParams() {

    // Vérification du type d'entrée
    if (!params.input_type) {
        error "ERREUR: Le paramètre --input_type est obligatoire ('bcl' ou 'fastq')."
    }
    if (!['bcl', 'fastq'].contains(params.input_type)) {
        error "ERREUR: --input_type doit être 'bcl' ou 'fastq'. Valeur reçue : '${params.input_type}'."
    }

    // Vérification du dossier d'entrée
    if (!params.input_dir) {
        error "ERREUR: Le paramètre --input_dir est obligatoire."
    }
    def input_path = file(params.input_dir)
    if (!input_path.exists()) {
        error "ERREUR: Le dossier d'entrée n'existe pas : ${params.input_dir}"
    }
    if (!input_path.isDirectory()) {
        error "ERREUR: --input_dir doit être un dossier, pas un fichier : ${params.input_dir}"
    }

    // Vérification du génome de référence
    if (!params.genome_reference) {
        error "ERREUR: Le paramètre --genome_reference est obligatoire."
    }
    def genome_path = file(params.genome_reference)
    if (!genome_path.exists()) {
        error "ERREUR: Le dossier de référence génomique n'existe pas : ${params.genome_reference}"
    }

    // Vérification de la sample sheet pour le mode BCL
    if (params.input_type == 'bcl') {
        if (!params.sample_sheet) {
            error "ERREUR: --sample_sheet est obligatoire en mode BCL."
        }
        def ss_path = file(params.sample_sheet)
        if (!ss_path.exists()) {
            error "ERREUR: La sample sheet n'existe pas : ${params.sample_sheet}"
        }
    }

    // Vérification du dossier de sortie
    if (!params.output_dir) {
        error "ERREUR: Le paramètre --output_dir est obligatoire."
    }

    log.info "✔ Validation des paramètres OK."
}

// ========================================================================================
// WORKFLOW PRINCIPAL
// ========================================================================================

workflow {

    printHeader()
    validateParams()

    // Channel pour les fichiers de version (pour MultiQC et traçabilité)
    ch_versions = Channel.empty()

    // -----------------------------------------------------------------------
    // CAS 1 : ENTRÉE BCL → mkfastq → count → multiqc
    // -----------------------------------------------------------------------
    if (params.input_type == 'bcl') {

        log.info "Mode : BCL détecté. Lancement de mkfastq puis cellranger count."

        // Channel : dossier BCL + sample sheet
        ch_bcl_input = Channel.of([
            file(params.input_dir, checkIfExists: true),
            file(params.sample_sheet, checkIfExists: true)
        ])

        // Sous-workflow BCL → COUNT
        BCL_TO_COUNT(ch_bcl_input)

        ch_versions = ch_versions.mix(BCL_TO_COUNT.out.versions)

        // Collecte des métriques pour MultiQC
        ch_qc_inputs = BCL_TO_COUNT.out.metrics
            .collect()
            .map { files -> files }

        MULTIQC(ch_qc_inputs)

    // -----------------------------------------------------------------------
    // CAS 2 : ENTRÉE FASTQ → count → multiqc
    // -----------------------------------------------------------------------
    } else if (params.input_type == 'fastq') {

        log.info "Mode : FASTQ détecté. Détection automatique des samples."

        // Détection automatique des FASTQ organisés par sample
        // Convention : <input_dir>/<sample_id>/*_R1_*.fastq.gz
        ch_fastq_dirs = Channel
            .fromPath("${params.input_dir}/*", type: 'dir')
            .map { dir ->
                def sample_id = dir.name
                def fastq_files = file("${dir}/*.fastq.gz")
                if (fastq_files.size() == 0) {
                    log.warn "ATTENTION: Aucun fichier FASTQ trouvé dans ${dir} — sample ignoré."
                    return null
                }
                return tuple(sample_id, dir)
            }
            .filter { it != null }

        // Vérification qu'au moins un sample a été trouvé
        ch_fastq_dirs
            .ifEmpty {
                error "ERREUR: Aucun dossier de sample FASTQ trouvé dans ${params.input_dir}. Vérifiez la structure des dossiers."
            }

        // Sous-workflow FASTQ → COUNT
        FASTQ_TO_COUNT(ch_fastq_dirs)

        ch_versions = ch_versions.mix(FASTQ_TO_COUNT.out.versions)

        // Collecte des métriques pour MultiQC
        ch_qc_inputs = FASTQ_TO_COUNT.out.metrics
            .collect()
            .map { files -> files }

        MULTIQC(ch_qc_inputs)
    }

    // Émission des versions pour reproductibilité
    ch_versions
        .unique()
        .collectFile(name: "${params.output_dir}/pipeline_versions.txt", newLine: true)

    log.info "✔ Pipeline terminé. Résultats dans : ${params.output_dir}"
}

// ========================================================================================
// GESTION DES ERREURS GLOBALES
// ========================================================================================

workflow.onError {
    log.error """
    ╔══════════════════════════════════════════════════════╗
    ║  PIPELINE ÉCHOUÉ — Consultez les logs pour détails  ║
    ╚══════════════════════════════════════════════════════╝
    Erreur : ${workflow.errorMessage}
    Rapport : ${workflow.launchDir}/pipeline_report.html
    """.stripIndent()
}

workflow.onComplete {
    if (workflow.success) {
        log.info """
    ╔══════════════════════════════════════════════════════╗
    ║           PIPELINE TERMINÉ AVEC SUCCÈS              ║
    ╚══════════════════════════════════════════════════════╝
    Durée totale   : ${workflow.duration}
    Résultats      : ${params.output_dir}
    Rapport HTML   : ${workflow.launchDir}/pipeline_report.html
    Timeline       : ${workflow.launchDir}/pipeline_timeline.html
    DAG            : ${workflow.launchDir}/pipeline_dag.html
        """.stripIndent()
    }
}
