# CellRanger_Flow

> Nextflow DSL2 pipeline for single-cell 10x Genomics processing on HPC clusters with SLURM scheduler.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%A9%A4-blue.svg)](https://www.nextflow.io/)
[![DSL2](https://img.shields.io/badge/DSL-2-brightgreen.svg)](https://www.nextflow.io/docs/latest/dsl2.html)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

---

## Overview

**CellRanger_Flow** is a production-ready Nextflow pipeline designed to automate the complete BCL-to-quantification workflow for 10x Genomics single-cell experiments on high-performance computing (HPC) clusters with SLURM scheduler support.

The pipeline orchestrates the following steps:

1. **Sample Sheet Preprocessing** – Validates and standardizes input sample sheet CSV format
2. **BCL-to-FASTQ Conversion** – Converts raw sequencing BCL files to FASTQ using `cellranger mkfastq`
3. **Multi-Modal Alignment** – Performs parallel batch-wise analysis with `cellranger multi` (Gene Expression + VDJ)
4. **QC Aggregation** – Generates unified quality control reports with MultiQC
5. **Version Tracking** – Consolidates tool versions from all steps for reproducibility

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Execution on HPC (SLURM)](#execution-on-hpc-slurm)
  - [Safe Testing with Stub Runs](#safe-testing-with-stub-runs)
  - [Long-Running Jobs](#long-running-jobs)
- [Input Parameters](#input-parameters)
- [Outputs](#outputs)
- [Project Structure](#project-structure)
- [Modules](#modules)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [Authors](#authors)

---

## Features

✅ **Automated BCL-to-Multi Workflow**
- End-to-end single-cell 10x processing pipeline
- Support for multi-modal experiments (GEX + VDJ)

✅ **HPC-Optimized**
- SLURM scheduler integration for cluster execution
- Per-batch parallel processing
- Configurable resource allocation

✅ **Robust Error Handling**
- Input validation at startup
- Verbose logging and error reporting
- Safe dry-run testing with stub blocks

✅ **Production-Ready**
- Nextflow best practices (DSL2, modular design)
- Comprehensive parameter validation
- Execution reports and traceability

---

## Prerequisites

### System Requirements

- **Nextflow** ≥ 23.10.0
- **Java** ≥ 17 (OpenJDK 21 recommended)
- **SLURM** scheduler (for HPC cluster submission)
- **Cell Ranger** 7.1.0+ (must be pre-installed on cluster)
- **MultiQC** 1.21+
- **Conda** (for environment management)

### Data Requirements

- Valid BCL directory from Illumina sequencer
- Sample sheet CSV with proper 10x format
- Pre-built Cell Ranger references (GEX and VDJ)

---

## Installation

### 1. Create Nextflow Conda Environment

```bash
# Load Conda module on HPC cluster
module load miniconda3

# Create dedicated environment
conda create -n nextflow_env \
  -c conda-forge -c bioconda \
  nextflow=25.10.4 openjdk=21 -y

# Activate environment
conda activate nextflow_env

# Verify
nextflow --version && java -version
```

### 2. Verify Prerequisites

```bash
# Cell Ranger
/labos/UGM/dev/cellranger-7.1.0/bin/cellranger --version

# MultiQC (in container or conda)
multiqc --version
```

---

## Quick Start

### Minimal Example

```bash
cd /home/dutel/cellranger-flow

# Load miniconda module
module load devel/python/Miniconda3-py39_4.10.3

# Source conda
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate environment
conda activate nextflow_env

# Test with stub-run (no computation, fast validation)
nextflow run main.nf \
  -c assets/nextflow.config \
  -params-file assets/params.json \
  -profile slurm \
  -stub-run

# Run actual pipeline
nextflow run main.nf \
  -c assets/nextflow.config \
  -params-file assets/params.json \
  -profile slurm \
  -resume
```

---

## Usage

### Execution on HPC (SLURM)

#### Option 1: SLURM Job Controller (Recommended)

Submit to SLURM for persistence even if you disconnect:

```bash
sbatch Slurm_job_nextflow.sh
```

Monitor:

```bash
squeue -u $USER
tail -f logs/nextflow_controller_*.out
```

#### Option 2: Persistent Terminal (tmux)

Keep job running via detachable terminal:

```bash
tmux new -s cellranger
nextflow run main.nf -c assets/nextflow.config -params-file assets/params.json -profile slurm

# Detach: Ctrl+b, d
# Re-attach: tmux attach -t cellranger
```

### Safe Testing with Stub Runs

Validate workflow logic without expensive computation:

```bash
nextflow run main.nf \
  -c assets/nextflow.config \
  -params-file assets/params.json \
  -profile slurm \
  -stub-run \
  -qs 1 \
  -process.maxForks 1
```

**Expected result:** Completes in seconds with stub files, no SLURM jobs submitted.

### Long-Running Jobs

**Key advantage:** SLURM keeps your Nextflow orchestrator running even if you disconnect.

```bash
# Submit and get job ID
sbatch Slurm_job_nextflow.sh  
# Output: Submitted batch job 12345678

# Later, check status (can reconnect later)
squeue -j 12345678
tail -f logs/nextflow_controller_12345678.out

# Cancel if needed
scancel 12345678
```

---

## Input Parameters

### Required Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `run_id` | String | Unique run identifier | `HCHNTDMX2` |
| `raw_sample_sheet_file_path` | Path | 10x sample sheet CSV | `/path/to/sample.csv` |
| `bcl_dir` | Path | BCL directory from sequencer | `/sequenceurs/NovaSeq1/...` |
| `batch_ids` | String | Comma-separated batch IDs | `74,75,76` |
| `protocol_prefix` | String | Batch name prefix | `MIDAS2` |
| `path_ref_gex` | Path | Cell Ranger GEX reference | `/path/to/refdata-gex-...` |
| `path_ref_vdj` | Path | Cell Ranger VDJ reference | `/path/to/refdata-vdj-...` |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `output_dir` | `${launchDir}/results/${run_id}` | Output base directory |
| `log_dir` | `${launchDir}/logs/${run_id}` | Logs base directory |
| `cpu_limit` | `16` | CPUs per process |
| `memory_limit` | `64` | Memory per process (GB) |
| `time_limit` | `168` | Time limit per process (hours) |
| `slurm_partition` | `phoenix` | SLURM partition name |

---

## Outputs

### Directory Structure

```
results/{run_id}/
├── preprocessed_sample_sheet/
│   └── Index_mkfastq_{run_id}.csv
├── qc/
│   └── fastq_output_{run_id}/
│       └── lane1/*.fastq.gz
├── alignment/
│   └── {batch_id}/outs/
│       ├── metrics_summary.csv
│       ├── web_summary.html
│       └── (Cell Ranger outputs)
└── multiqc/
    ├── multiqc_report_{run_id}.html
    └── multiqc_data_{run_id}/

logs/{run_id}/
├── preprocessing/*.log
├── qc/*.log
├── alignment/*.log
├── multiqc/*.log
├── {date}_versions.yaml
└── {date}_nextflow.log
```

### Key Outputs

| File | Purpose | Location |
|------|---------|----------|
| `Index_mkfastq_*.csv` | Standardized sample sheet | `preprocessing_output_dir` |
| `fastq_output_*/` | FASTQ files (per lane) | `qc_output_dir` |
| `*/outs/metrics_summary.csv` | QC metrics per batch | `alignment_output_dir` |
| `multiqc_report_*.html` | Unified QC report | `multiqc_output_dir` |
| `*_versions.yaml` | Tool versions for reproducibility | `log_dir` |

---

## Project Structure

```
cellranger-flow/
├── main.nf                      # Main workflow orchestration (DSL2)
├── README.md                    # Documentation (this file)
├── Slurm_job_nextflow.sh        # SLURM submission script
│
├── assets/
│   ├── nextflow.config          # Global configuration
│   ├── params.json              # Parameter file
│   ├── logging.sh               # Logging functions
│   └── dag-Entry.html           # Example DAG
│
├── modules/                     # Reusable DSL2 modules
│   ├── preprocessing.nf         # Sample sheet normalization
│   ├── cellranger_mkfastq.nf    # BCL → FASTQ conversion
│   ├── cellranger_multi.nf      # Alignment (parallel per batch)
│   ├── multiqc.nf               # QC aggregation
│   └── merge_versions.nf        # Version consolidation
│
└── subworkflows/                # Future workflows
```

---

## Modules

| Module | Task | Parallelization |
|--------|------|-----------------|
| **PREPROCESSING** | Standardize sample sheet | Single |
| **CELLRANGER_MKFASTQ** | Convert BCL to FASTQ | Single |
| **CELLRANGER_MULTI** | Alignment (GEX + VDJ) | **Parallel (per batch)** |
| **MULTIQC** | Aggregate QC reports | Single |
| **MERGE_VERSIONS** | Consolidate versions | Single |

---

## Troubleshooting

### Issue: `-stub-run` still submits SLURM jobs

**Cause:** Stub blocks not defined in modules.

**Solution:** Verify all modules define `stub: { }`:

```bash
grep "stub:" modules/*.nf
```

### Issue: File name collision on `versions.yml`

**Cause:** Multiple modules emit identical filename.

**Solution:** Use unique per-module names:
- `1_preprocessing_versions.yml`
- `2_cellranger_mkfastq_versions.yml`
- `3_cellranger_multi_versions.yml`
- `4_multiqc_versions.yml`

### Issue: "cellranger: command not found"

**Solution:**

```bash
# Add to PATH
export PATH="/labos/UGM/dev/cellranger-7.1.0/bin:$PATH"
cellranger --version

# Or update Slurm script with module load
```

### Issue: "Java version requirement not met"

**Solution:** Ensure Java 17+ in Nextflow environment:

```bash
conda create -n nextflow_env \
  openjdk=21 nextflow=25.10.4 \
  -c conda-forge
```

### Issue: SLURM partition doesn't exist

**Check available partitions:**

```bash
sinfo -a
```

**Update config:**

```bash
params.slurm_partition = "your_partition"
```

---

## Citation

If you use CellRanger_Flow in your research, please cite:

```bibtex
@software{cellranger_flow_2026,
  title={CellRanger_Flow: Nextflow DSL2 for 10x Genomics Single-Cell Processing},
  author={Dutel, Jordan and GENIM Team},
  year={2026},
  url={https://github.com/your-org/cellranger-flow}
}
```

---

## Authors

- **Jordan Dutel** – GENIM Team, CRCT, INSERM
  - Email: jordan.dutel@inserm.fr
  - Role: Pipeline Developer & Maintainer

### Contributing

Contributions welcome! Please:
1. Create a feature branch
2. Commit with clear messages
3. Submit a pull request
4. Ensure all validations pass

---

## License

MIT License – see [LICENSE](LICENSE) file for details.

---

## Support

For issues, questions, or feedback:

1. Check [Troubleshooting](#troubleshooting) section above
2. Review pipeline logs (`.nextflow.log`, `logs/nf_controller_*.out`)
3. Verify parameters in `assets/params.json`
4. Open an issue on GitHub with:
   - Command used
   - Error messages
   - `.nextflow.log` excerpt
   - System/cluster info

---

**Last Updated:** April 2026  
**Pipeline Version:** 1.0.0  
**Nextflow Version:** ≥25.10.4  
**Cell Ranger Version:** ≥7.1.0  
**Maintained by:** GENIM Team, CRCT

