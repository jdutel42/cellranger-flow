# sc10x_pipeline — Single-Cell 10x Genomics Pipeline (DSL2)

Pipeline Nextflow DSL2 pour le traitement automatisé de données **Single-Cell 10x Genomics**,
depuis les fichiers BCL ou FASTQ jusqu'aux matrices de comptage et rapports QC.

---

## Table des matières

1. [Architecture du projet](#architecture)
2. [Prérequis](#prérequis)
3. [Installation](#installation)
4. [Paramètres](#paramètres)
5. [Exemples de commandes](#exemples)
6. [Structure des sorties](#sorties)
7. [Compatibilité nf-core](#nf-core)

---

## Architecture du projet <a name="architecture"></a>

```
sc10x_pipeline/
├── main.nf                          ← Point d'entrée principal
├── nextflow.config                  ← Configuration globale + profils
│
├── modules/                         ← Modules atomiques (1 process = 1 outil)
│   ├── cellranger_mkfastq.nf        ← BCL → FASTQ
│   ├── cellranger_count.nf          ← FASTQ → matrices + métriques
│   └── multiqc.nf                   ← Agrégation QC
│
├── subworkflows/                    ← Sous-workflows réutilisables
│   ├── bcl_to_count.nf              ← mkfastq → count (mode BCL)
│   └── fastq_to_count.nf            ← count direct (mode FASTQ)
│
├── assets/
│   ├── sample_sheet.csv             ← Exemple de sample sheet Cell Ranger
│   └── multiqc_config.yaml          ← Configuration MultiQC
│
└── README.md
```

### Flux de données

```
Mode BCL :
  [BCL dir] + [sample_sheet.csv]
      │
      ▼
  CELLRANGER_MKFASTQ
      │
      ▼ (un dossier FASTQ par sample)
  CELLRANGER_COUNT ×N  (parallèle)
      │
      ▼
  MULTIQC  (agrégation)

Mode FASTQ :
  [input_dir/sample_1/] [input_dir/sample_2/] ...
      │
      ▼ (détection automatique)
  CELLRANGER_COUNT ×N  (parallèle)
      │
      ▼
  MULTIQC  (agrégation)
```

---

## Prérequis <a name="prérequis"></a>

| Outil       | Version minimale | Rôle                        |
|-------------|------------------|-----------------------------|
| Nextflow    | 23.10.0          | Orchestrateur de pipeline   |
| Docker      | 20.x+            | Conteneurisation (local)    |
| Singularity | 3.8+             | Conteneurisation (HPC)      |
| Java        | 11+              | Runtime Nextflow            |

Les outils bioinformatiques (Cell Ranger, MultiQC) sont **embarqués dans les containers**
et ne nécessitent pas d'installation manuelle.

---

## Installation <a name="installation"></a>

```bash
# Cloner le dépôt
git clone https://github.com/myorg/sc10x_pipeline.git
cd sc10x_pipeline

# Vérifier la version de Nextflow
nextflow -version   # doit être >= 23.10.0

# Mettre à jour Nextflow si nécessaire
nextflow self-update
```

---

## Paramètres <a name="paramètres"></a>

| Paramètre          | Défaut                  | Description                                       |
|--------------------|-------------------------|---------------------------------------------------|
| `--input_type`     | `fastq`                 | `bcl` ou `fastq`                                  |
| `--input_dir`      | (obligatoire)           | Dossier BCL ou racine des FASTQ                   |
| `--output_dir`     | `./results`             | Dossier de sortie                                 |
| `--genome_reference` | (obligatoire)         | Chemin vers la référence Cell Ranger pré-buildée  |
| `--sample_sheet`   | (requis si BCL)         | Sample sheet CSV (format Cell Ranger)             |
| `--localcores`     | `16`                    | CPUs alloués à Cell Ranger                        |
| `--localmemory`    | `64`                    | RAM allouée en GB                                 |
| `--chemistry`      | `auto`                  | Chimie 10x (ex: `SC3Pv3`, `SC5P-PE`)             |
| `--expect_cells`   | `5000`                  | Nombre de cellules attendues par sample           |
| `--include_introns`| `true`                  | Inclure les reads introniques                     |
| `--force_cells`    | (désactivé)             | Forcer un nombre précis de cellules               |

---

## Exemples de commandes <a name="exemples"></a>

### Mode FASTQ — Docker (local)

```bash
nextflow run main.nf \
    -profile docker \
    --input_type fastq \
    --input_dir /data/fastq_samples \
    --genome_reference /references/refdata-gex-GRCh38-2020-A \
    --output_dir /results/run_2024_01 \
    --localcores 16 \
    --localmemory 64 \
    --expect_cells 8000 \
    --include_introns true
```

### Mode BCL — Singularity (local HPC)

```bash
nextflow run main.nf \
    -profile singularity \
    --input_type bcl \
    --input_dir /data/bcl/240115_A00123_0001_BHXXXXXX \
    --sample_sheet /data/samplesheets/run_240115.csv \
    --genome_reference /references/refdata-gex-GRCh38-2020-A \
    --output_dir /results/run_240115 \
    --localcores 32 \
    --localmemory 128
```

### Mode BCL — SLURM + Singularity (cluster HPC)

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --input_type bcl \
    --input_dir /scratch/user/bcl/240115_A00123 \
    --sample_sheet /scratch/user/sheets/run_240115.csv \
    --genome_reference /shared/references/refdata-gex-GRCh38-2020-A \
    --output_dir /scratch/user/results/run_240115 \
    --localcores 32 \
    --localmemory 256 \
    -resume
```

### Reprise d'un run interrompu

```bash
# L'option -resume permet de reprendre sans recalculer les étapes déjà terminées
nextflow run main.nf -profile docker --input_type fastq ... -resume
```

### Mode stub (test sans données réelles)

```bash
nextflow run main.nf \
    -profile docker \
    -stub \
    --input_type fastq \
    --input_dir /tmp/test_fastq \
    --genome_reference /tmp/test_genome \
    --output_dir /tmp/test_results
```

---

## Structure des sorties <a name="sorties"></a>

```
results/
├── cellranger_count/
│   ├── Sample_A/
│   │   └── outs/
│   │       ├── filtered_feature_bc_matrix/   ← Matrices (MEX format)
│   │       │   ├── matrix.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   └── features.tsv.gz
│   │       ├── metrics_summary.csv           ← Métriques QC
│   │       ├── web_summary.html              ← Rapport HTML Cell Ranger
│   │       └── molecule_info.h5             ← Pour cellranger aggr
│   └── Sample_B/
│       └── outs/
│           └── ...
│
├── mkfastq/                                  ← (mode BCL uniquement)
│   └── run_id_mkfastq/
│       └── fastq_output/
│           ├── Sample_A/
│           └── Sample_B/
│
├── multiqc/
│   ├── multiqc_report.html                  ← Rapport QC global
│   └── multiqc_data/                        ← Données brutes MultiQC
│
├── logs/
│   ├── mkfastq/
│   │   └── mkfastq_*.log
│   └── cellranger_count/
│       └── count_*.log
│
└── pipeline_versions.txt                    ← Versions de tous les outils
```

---

## Compatibilité nf-core <a name="nf-core"></a>

Ce pipeline suit les conventions nf-core pour faciliter son intégration :

### Structure compatible nf-core

```
sc10x_pipeline/
├── bin/                             ← Scripts helper (check_samplesheet.py, etc.)
├── conf/
│   ├── base.config                  ← Ressources par défaut
│   ├── igenomes.config              ← Références IGENOMES
│   └── modules.config               ← Options par module (ext.args)
├── docs/
│   ├── usage.md
│   └── output.md
├── lib/
│   ├── NfcoreSchema.groovy          ← Validation JSON schema
│   └── WorkflowSc10x.groovy        ← Fonctions métier
├── modules/
│   └── nf-core/                    ← Modules nf-core (si utilisés)
├── workflows/
│   └── sc10x.nf                    ← Workflow principal (convention nf-core)
├── assets/
│   ├── schema_input.json            ← Validation de la sample sheet
│   └── multiqc_config.yaml
├── main.nf
├── nextflow.config
├── nextflow_schema.json            ← Schéma de validation des paramètres
└── CITATIONS.md                    ← Citations des outils
```

### Adaptations recommandées pour nf-core

1. **`nextflow_schema.json`** : Valider les paramètres avec `nf-validation`.
2. **`conf/modules.config`** : Centraliser les `ext.args` de chaque module.
3. **`lib/WorkflowSc10x.groovy`** : Déplacer `validateParams()` et `printHeader()`.
4. **`bin/check_samplesheet.py`** : Script Python de validation de la sample sheet.
5. **`CITATIONS.md`** : Référencer Cell Ranger, MultiQC, Nextflow.

### Commande de validation nf-core

```bash
# Installer nf-core tools
pip install nf-core

# Valider la structure du pipeline
nf-core lint .

# Créer un pipeline à partir du template nf-core
nf-core create --name sc10x --description "Single-Cell 10x Genomics Pipeline"
```

---

## Outils et versions

| Outil       | Version | Container                                          |
|-------------|---------|-----------------------------------------------------|
| Cell Ranger | 7.2.0   | `nfcore/cellranger:7.2.0`                          |
| MultiQC     | 1.21    | `quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0` |
| Nextflow    | 23.10+  | (moteur local)                                      |

---

## Licence

MIT © 2024 Bioinformatics Pipeline Team
