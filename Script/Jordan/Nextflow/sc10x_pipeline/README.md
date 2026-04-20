# sc10x_pipeline - Single-Cell 10x Genomics Pipeline (DSL2)

Nextflow DSL2 pipeline for automated processing of **Single-Cell 10x Genomics** data,
from BCL or FASTQ files to count matrices and QC reports.

---

## Table of Contents

1. [Project Architecture](#architecture)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Parameters](#parameters)
5. [Command Examples](#examples)
6. [Output Structure](#outputs)
7. [nf-core Compatibility](#nf-core)

---

## Project Architecture <a name="architecture"></a>

```
sc10x_pipeline/
├── main.nf                          <- Main entry point
├── nextflow.config                  <- Global configuration + profiles
│
├── modules/                         <- Atomic modules (1 process = 1 tool)
│   ├── cellranger_mkfastq.nf        <- BCL -> FASTQ
│   ├── cellranger_count.nf          <- FASTQ -> matrices + metrics
│   └── multiqc.nf                   <- QC aggregation
│
├── subworkflows/                    <- Reusable sub-workflows
│   ├── bcl_to_count.nf              <- mkfastq -> count (BCL mode)
│   └── fastq_to_count.nf            <- direct count (FASTQ mode)
│
├── assets/
│   ├── sample_sheet.csv             <- Example Cell Ranger sample sheet
│   └── multiqc_config.yaml          <- MultiQC configuration
│
└── README.md
```

### Data Flow

```
BCL mode:
  [BCL dir] + [sample_sheet.csv]
      |
      v
  CELLRANGER_MKFASTQ
      |
      v (one FASTQ directory per sample)
  CELLRANGER_COUNT xN  (parallel)
      |
      v
  MULTIQC  (aggregation)

FASTQ mode:
  [input_dir/sample_1/] [input_dir/sample_2/] ...
      |
      v (automatic detection)
  CELLRANGER_COUNT xN  (parallel)
      |
      v
  MULTIQC  (aggregation)
```

---

## Prerequisites <a name="prerequisites"></a>

| Tool        | Minimum version | Role                         |
|-------------|------------------|------------------------------|
| Nextflow    | 23.10.0          | Pipeline orchestrator        |
| Docker      | 20.x+            | Containerization (local)     |
| Singularity | 3.8+             | Containerization (HPC)       |
| Java        | 11+              | Nextflow runtime             |

Bioinformatics tools (Cell Ranger, MultiQC) are **bundled in containers**
and do not require manual installation.

---

## Installation <a name="installation"></a>

```bash
# Clone repository
git clone https://github.com/myorg/sc10x_pipeline.git
cd sc10x_pipeline

# Check Nextflow version
nextflow -version   # must be >= 23.10.0

# Update Nextflow if needed
nextflow self-update
```

---

## Parameters <a name="parameters"></a>

| Parameter            | Default               | Description                                  |
|----------------------|-----------------------|----------------------------------------------|
| `--input_type`       | `fastq`               | `bcl` or `fastq`                             |
| `--input_dir`        | (required)            | BCL directory or FASTQ root directory        |
| `--output_dir`       | `./results`           | Output directory                              |
| `--genome_reference` | (required)            | Path to pre-built Cell Ranger reference      |
| `--sample_sheet`     | (required if BCL)     | CSV sample sheet (Cell Ranger format)        |
| `--localcores`       | `16`                  | CPUs allocated to Cell Ranger                |
| `--localmemory`      | `64`                  | RAM allocated in GB                          |
| `--chemistry`        | `auto`                | 10x chemistry (e.g. `SC3Pv3`, `SC5P-PE`)     |
| `--expect_cells`     | `5000`                | Expected number of cells per sample          |
| `--include_introns`  | `true`                | Include intronic reads                       |
| `--force_cells`      | (disabled)            | Force a specific number of cells             |

---

## Command Examples <a name="examples"></a>

### FASTQ mode - Docker (local)

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

### BCL mode - Singularity (local HPC)

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

### BCL mode - SLURM + Singularity (HPC cluster)

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

### Resume an interrupted run

```bash
# The -resume option continues without recomputing completed steps
nextflow run main.nf -profile docker --input_type fastq ... -resume
```

### Stub mode (test without real data)

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

## Output Structure <a name="outputs"></a>

```
results/
├── cellranger_count/
│   ├── Sample_A/
│   │   └── outs/
│   │       ├── filtered_feature_bc_matrix/   <- Matrices (MEX format)
│   │       │   ├── matrix.mtx.gz
│   │       │   ├── barcodes.tsv.gz
│   │       │   └── features.tsv.gz
│   │       ├── metrics_summary.csv           <- QC metrics
│   │       ├── web_summary.html              <- Cell Ranger HTML report
│   │       └── molecule_info.h5              <- For cellranger aggr
│   └── Sample_B/
│       └── outs/
│           └── ...
│
├── mkfastq/                                  <- (BCL mode only)
│   └── run_id_mkfastq/
│       └── fastq_output/
│           ├── Sample_A/
│           └── Sample_B/
│
├── multiqc/
│   ├── multiqc_report.html                   <- Global QC report
│   └── multiqc_data/                         <- Raw MultiQC data
│
├── logs/
│   ├── mkfastq/
│   │   └── mkfastq_*.log
│   └── cellranger_count/
│       └── count_*.log
│
└── pipeline_versions.txt                     <- Versions of all tools
```

---

## nf-core Compatibility <a name="nf-core"></a>

This pipeline follows nf-core conventions to facilitate integration.

### nf-core Compatible Structure

```
sc10x_pipeline/
├── bin/                             <- Helper scripts (check_samplesheet.py, etc.)
├── conf/
│   ├── base.config                  <- Default resources
│   ├── igenomes.config              <- IGENOMES references
│   └── modules.config               <- Module-level options (ext.args)
├── docs/
│   ├── usage.md
│   └── output.md
├── lib/
│   ├── NfcoreSchema.groovy          <- JSON schema validation
│   └── WorkflowSc10x.groovy         <- Workflow helper functions
├── modules/
│   └── nf-core/                     <- nf-core modules (if used)
├── workflows/
│   └── sc10x.nf                     <- Main workflow (nf-core convention)
├── assets/
│   ├── schema_input.json            <- Sample sheet validation
│   └── multiqc_config.yaml
├── main.nf
├── nextflow.config
├── nextflow_schema.json             <- Parameter validation schema
└── CITATIONS.md                     <- Tool citations
```

### Recommended nf-core Adaptations

1. **`nextflow_schema.json`**: Validate parameters with `nf-validation`.
2. **`conf/modules.config`**: Centralize `ext.args` for each module.
3. **`lib/WorkflowSc10x.groovy`**: Move `validateParams()` and `printHeader()` there.
4. **`bin/check_samplesheet.py`**: Add Python sample sheet validation.
5. **`CITATIONS.md`**: Reference Cell Ranger, MultiQC, and Nextflow.

### nf-core Validation Command

```bash
# Install nf-core tools
pip install nf-core

# Validate pipeline structure
nf-core lint .

# Create a pipeline from nf-core template
nf-core create --name sc10x --description "Single-Cell 10x Genomics Pipeline"
```

---

## Tools and Versions

| Tool        | Version | Container                                          |
|-------------|---------|----------------------------------------------------|
| Cell Ranger | 7.2.0   | `nfcore/cellranger:7.2.0`                          |
| MultiQC     | 1.21    | `quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0` |
| Nextflow    | 23.10+  | (local engine)                                     |

---

## License

MIT (c) 2024 Bioinformatics Pipeline Team
