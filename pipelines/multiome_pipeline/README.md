# Multiome Pipeline

A comprehensive pipeline for analyzing 10x Multiome (RNA + ATAC) single-cell data. This pipeline processes raw multiome data through quality control, clustering, cell type annotation, and exports results in H5 format for visualization in ISCVAM.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Quick Start](#quick-start)
- [Pipeline Steps](#pipeline-steps)
- [Configuration](#configuration)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)

## Prerequisites

- **R** >= 4.0
- **Conda** (Miniconda or Anaconda) - for MACS3 peak calling
- **10x Multiome data** in standard Cell Ranger ARC output format

## Installation

### 1. Install R Packages and MACS3

Run the installation script from the `multiome_pipeline` directory:

```r
setwd("pipelines/multiome_pipeline")
source("code_organized/install_packages.R")
```

This script will:
- Install required R packages (Seurat, Signac, SingleR, etc.)
- Create a conda environment with MACS3 for peak calling

### 2. Required R Packages

| Package | Purpose |
|---------|---------|
| Seurat | Single-cell analysis |
| Signac | ATAC-seq analysis |
| SingleR | Automated cell type annotation |
| EnsDb.Hsapiens.v86 | Human gene annotations |
| BSgenome.Hsapiens.UCSC.hg38 | Human genome sequence |
| dplyr | Data manipulation |
| future | Parallel processing |
| rhdf5 | HDF5 file handling |

### 3. Reference Data

Ensure the SingleR references are available:
```
code_organized/ref_data/singler.refs.rdata
```

## Directory Structure

```
multiome_pipeline/
├── code_organized/
│   ├── configure_dataset.R      # Dataset configuration (edit this)
│   ├── install_packages.R       # Installation script
│   ├── core/
│   │   ├── run_pipeline.R                # Main pipeline script
│   │   └── multiome_functions_pipeline.R # Core analysis functions
│   ├── utils/
│   │   ├── iscvamR_h5_customized.SeuratV5.R      # H5 utilities
│   │   └── updated_write_h5_scripts_using_iscvamR.R # H5 export
│   └── ref_data/
│       └── singler.refs.rdata            # SingleR references
├── example_data/
│   └── <dataset_name>/
│       └── raw_data/            # Place your 10x data here
└── results/
    └── results_<dataset_name>/  # Output files
```

## Quick Start

### Step 1: Prepare Your Data

Place your 10x Multiome data in the `example_data` folder:
```
example_data/<your_dataset_name>/raw_data/
├── barcodes.tsv.gz
├── features.tsv.gz
├── matrix.mtx.gz
├── atac_fragments.tsv.gz
├── atac_fragments.tsv.gz.tbi
└── atac_peak_annotation.tsv
```

### Step 2: Configure the Dataset

Edit `code_organized/configure_dataset.R`:

```r
# Set your dataset name and species
dataset_name <- "your_dataset_name"
species <- "Human"  # or "Mouse"
```

### Step 3: Run the Pipeline

From R, in the `multiome_pipeline` directory:

```r
# Set working directory
setwd("pipelines/multiome_pipeline")

# Load configuration
source("code_organized/configure_dataset.R")

# Run the pipeline
source("code_organized/core/run_pipeline.R")
```

## Pipeline Steps

The pipeline executes the following steps:

| Step | Description | Output |
|------|-------------|--------|
| 1 | Create Seurat object from raw data | `seurat.raw.RData` |
| 2 | Set up results directory | - |
| 3 | SingleR cell type annotation | `singler.RData` |
| 4 | Peak calling with MACS3 | `seurat.peaks.RData` |
| 5 | Quality control (RNA & ATAC) | `seurat.qc_rna.RData`, `seurat.qc_atac.RData` |
| 6 | RNA analysis & clustering | `seurat.analyzed.rna.RData`, `clustering.rna.RData` |
| 7 | ATAC analysis & clustering | `seurat.analyzed.atac.RData`, `clustering.atac.RData` |
| 8 | WNN integration | `seurat.wnn.RData`, `seurat.wnn.clusterings.RData` |
| 8.1 | Azimuth annotation (optional) | Updated clusterings |
| 9 | Find markers | `markers_*_clusterings.RData` |
| 10 | Cell type statistics | `cell.and.cluster.stats.RData`, `covs.RData` |
| 11 | ISCVAM artifacts | `heatmap_artifacts.RData` |
| 12 | Export H5 file | `<dataset_name>_analysis.h5` |

## Configuration

### Resolution Settings

Clustering resolutions are defined in the functions file. Default resolutions:
```r
resolution.lst <- c(0.2, 0.4, 0.6, 0.8, 1.0)
```

### Azimuth References

Available Azimuth references for cell type annotation:
- `pbmcref` - PBMC
- `kidneyref` - Human Kidney
- `humancortexref` - Human Brain Cortex

Install references before use:
```r
SeuratData::InstallData("kidneyref")
```

```

## Output Files

### Results Directory

All output files are saved to `results/results_<dataset_name>/`:

| File | Description |
|------|-------------|
| `*.RData` | Intermediate R objects for each step |
| `cells.and.clusterings.xlsx` | Cell and cluster statistics |
| `<dataset_name>_analysis.h5` | Final H5 file for ISCVAM |

### Final H5 File

The H5 file contains:
- Normalized expression matrices (RNA & ATAC)
- UMAP coordinates (RNA, ATAC, WNN)
- Clustering assignments at multiple resolutions
- Cell type annotations (SingleR, Azimuth)
- Marker genes/peaks for each cluster
- Heatmap data for visualization

## Troubleshooting

### MACS3 Not Found

If peak calling fails:
```r
# Check if MACS3 is installed
system("conda run -n macs3 macs3 --version")

# Reinstall if needed
system("conda create -n macs3 python=3.9 -y")
system("conda run -n macs3 pip install MACS3")
```

### Memory Issues

For large datasets, increase R memory limits:
```r
options(future.globals.maxSize = 16 * 1024^3)  # 16GB
```

### Missing SingleR References

Download and save references:
```r
# Create references
library(celldex)
refs <- list(
  hpca = HumanPrimaryCellAtlasData(),
  blueprint = BlueprintEncodeData()
)
save(refs, file = "code_organized/ref_data/singler.refs.rdata")
```

### Dataset Not Found

Ensure your data directory structure matches:
```
example_data/<dataset_name>/raw_data/
```

And update `configure_dataset.R`:
```r
dataset_name <- "<your_dataset_name>"
```

## Support

For issues or questions, please refer to the main ISCVAM documentation or contact the development team.
