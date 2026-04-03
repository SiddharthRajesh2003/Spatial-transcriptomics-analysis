# Spatial Transcriptomics Analysis Pipeline — Glioblastoma

End-to-end pipeline for processing and analyzing spatial and single-cell RNA-seq data from glioblastoma (GBM) patient samples. Integrates a Nextflow preprocessing pipeline with a downstream Python analysis workflow using Cell2location, Squidpy, and scVI.

---

## Project Overview

This pipeline processes data from **3 primary** and **1 recurrent** GBM patient samples (GSE214966, scRNA-seq) alongside **1 Visium CytAssist** spatial transcriptomics slide ([10x Genomics Dataset](https://www.10xgenomics.com/datasets/gene-and-protein-expression-library-of-human-glioblastoma-cytassist-ffpe-2-standard)). The goal is to map cell-cell interactions in the GBM tumor microenvironment at spatial resolution.

**Key analyses:**
- Cell type annotation (Tumor, Microglia, TAM subtypes, CD8 T cells, OPC, Oligodendrocytes, Endothelial)
- Spatial deconvolution via Cell2location
- Neighborhood enrichment and co-occurrence analysis
- Ligand-receptor interaction analysis (Squidpy)
- Spatially variable gene detection (Moran's I)

---

## Repository Structure

```
.
├── main.nf                  # Nextflow pipeline entry point
├── nextflow.config          # Pipeline parameters and SLURM configuration
├── run_pipeline.sh          # SBATCH submission script
├── env.yaml                 # Conda environment for downstream analysis
├── code.qmd                 # Main downstream analysis (scRNA-seq + spatial)
├── downstream_analysis.qmd  # Additional downstream analysis
├── modules/
│   ├── fastqc/              # FastQC module
│   ├── cellranger/          # Cell Ranger module (scRNA-seq)
│   ├── multiqc/             # MultiQC module
│   └── spaceranger/         # Space Ranger module (Visium)
├── metadata/                # Sample metadata files
├── misc/                    # Helper scripts (samplesheet generation)
├── h5ad_files/              # Saved AnnData checkpoints
└── results/                 # Pipeline outputs
    ├── cellranger/          # Cell Ranger outputs per sample
    ├── spaceranger/         # Space Ranger outputs
    ├── cellbender/          # CellBender ambient RNA corrected outputs
    ├── fastqc/              # FastQC reports
    └── multiqc/             # Aggregated QC report
```

---

## Pipeline

### Stage 1 — Nextflow Preprocessing (`main.nf`)

Runs on HPC (SLURM) using Apptainer containers.

**Steps:**
1. **FastQC** — Quality control on raw FASTQ files for all samples
2. **Cell Ranger** — Alignment and UMI counting for scRNA-seq samples (GRCh38-2024-A reference)
3. **Space Ranger** — Alignment and spatial barcode detection for Visium CytAssist slide
4. **MultiQC** — Aggregated QC report across all samples

**Run the pipeline:**
```bash
sbatch run_pipeline.sh
```

**Skip completed steps** (e.g. if Cell Ranger already ran):
```bash
nextflow run main.nf -profile slurm --skip_cellranger --skip_fastqc
```

**Samplesheet formats:**

`cellranger_samples.csv`:
```
sample_id,fastq_dir
primary_001,/path/to/fastqs/primary_001
```

`spaceranger_samples.csv`:
```
fastqs,sample,library_type
/path/to/fastqs,GBM_CytAssist,Gene Expression
```

---

### Stage 2 — Downstream Analysis (`code.qmd`)

Run locally or on HPC after Nextflow completes. Reads CellBender-corrected outputs from `results/cellbender/`.

**Steps:**
1. Load and QC scRNA-seq data (CellBender + Scrublet doublet removal)
2. Batch integration with scVI
3. Leiden clustering + cell type annotation
4. TAM sub-clustering (TAM_APOE, TAM_SPP1, TAM_macrophage, Monocyte)
5. Load and QC Visium spatial data
6. Cell2location RegressionModel (cell-type signatures from scRNA-seq)
7. Cell2location spatial deconvolution
8. Neighborhood enrichment + co-occurrence analysis
9. Ligand-receptor analysis (Squidpy, 1000 permutations)
10. Spatial autocorrelation — Moran's I

**Render the analysis:**
```bash
quarto render code.qmd
```

---

## Environment Setup

```bash
conda env create -f env.yaml
conda activate scrna_env
```

**Key dependencies:**
- `scanpy` — single-cell analysis
- `scvi-tools` — batch integration
- `cell2location` — spatial deconvolution
- `squidpy` — spatial analysis and LR testing
- `scrublet` — doublet detection
- `rapids-singlecell` — GPU-accelerated single-cell analysis (requires CUDA 12.8)

---

## Checkpointing

Intermediate results are saved to `h5ad_files/` to avoid recomputing expensive steps:

| File | Contents |
|---|---|
| `scvi_ref.h5ad` | scVI-integrated scRNA-seq object |
| `annotated.h5ad` | Cell-type annotated scRNA-seq object |
| `ref_regression_c2l.h5ad` | Cell2location RegressionModel output |
| `c2l_complete.h5ad` | Cell2location spatial deconvolution output |
| `sp_protein.h5ad` | Antibody capture (protein) spatial object |

---

## HPC Configuration

| Process | CPUs | Memory | Time |
|---|---|---|---|
| FastQC | 4 | 8 GB | 1 h |
| Cell Ranger | 16 | 75 GB | 24 h |
| Space Ranger | 16 | 75 GB | 24 h |
| MultiQC | 4 | 16 GB | 6 h |
| Pipeline (head node) | 12 | 10 GB | 48 h |

---

## Singularity / Apptainer

See [singularity.md](singularity.md) for instructions on building and using the Apptainer containers for each module.

---

## Author

Siddharth Rajesh — Krolab, Indiana University
