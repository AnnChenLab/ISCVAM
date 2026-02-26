# scRNA-seq Pipeline

A comprehensive pipeline for analyzing single-cell RNA sequencing data and exporting results in H5 format for visualization in ISCVAM.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Quick Start](#quick-start)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Customization](#customization)

## Overview

This pipeline processes scRNA-seq data through two stages:

1. **First Stage (All Cells)**: QC, normalization, dimensionality reduction, clustering at multiple resolutions, marker discovery, and cell type annotation
2. **Second Stage (Per-Layer)**: Re-analyzes cell subsets (e.g., by major cell type) at higher resolution for refined structure

The final output is an ISCVAM-compatible H5 file with multi-layer analysis results.

## Prerequisites

- **R** >= 4.0
- **TISCH data** (or compatible 10x Genomics scRNA-seq format)

### Required R Packages

| Package | Purpose |
|---------|---------|
| Seurat | Single-cell analysis |
| SingleR | Automated cell type annotation |
| rhdf5 | HDF5 file handling |
| dplyr | Data manipulation |
| future | Parallel processing |
| here | Project-relative paths |
| tibble | Data frames |
| janitor | Data cleaning |

## Installation

Install required R packages:

```r
# CRAN packages
install.packages(c("dplyr", "future", "here", "tibble", "janitor"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Seurat", "SingleR", "rhdf5"))
```

## Directory Structure

```
sc_rna_pipeline/
├── README.md                   # This file
├── code/
│   ├── pipeline.R              # Main pipeline script
│   ├── scRNA_functions_pipeline.R  # Core analysis functions
│   ├── iscvamR_h5_customized.R     # H5 export utilities
│   └── pipeline.md             # Detailed pipeline documentation
├── reference_data/
│   ├── singler.refs.rdata      # SingleR reference datasets
│   └── meta.data_check.endo.track.progress.man.csv  # Dataset metadata
└── example.results/            # Output directory
    └── <dataset_name>/         # Per-dataset results
```

## Data and Metadata Setup

### Data Directory Structure

The pipeline expects TISCH-formatted scRNA-seq data. Set `data.path` in `pipeline.R` to point to your data root directory:

```r
data.path <- "/path/to/your/TISCH2_download"
```

Each dataset should be organized as a subdirectory containing:

```
<data.path>/
├── <Dataset.Name>/
│   ├── *.h5                           # 10x Genomics count matrix (HDF5)
│   └── *_CellMetainfo_table.tsv       # TISCH cell metadata table
├── <Another_Dataset>/
│   ├── *.h5
│   └── *_CellMetainfo_table.tsv
└── ...
```

**Example:**
```
/path/to/TISCH2_download/
├── NSCLC_GSE127471/
│   ├── NSCLC_GSE127471_expression.h5
│   └── NSCLC_GSE127471_CellMetainfo_table.tsv
├── Glioma_GSE84465/
│   ├── Glioma_GSE84465_expression.h5
│   └── Glioma_GSE84465_CellMetainfo_table.tsv
└── ...
```

### Dataset Metadata CSV

The pipeline uses a metadata CSV file to manage multiple datasets. This file is located at:

```
reference_data/meta.data_check.endo.track.progress.man.csv
```

Load it in your pipeline:
```r
meta.tisch <- read.csv(here("reference_data/meta.data_check.endo.track.progress.man.csv"))
```

#### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `Dataset.Name` | Unique dataset identifier (must match folder name in data.path) | `NSCLC_GSE127471` |
| `ID` | TISCH internal ID | `T010059` |
| `Species` | Human or Mouse | `Human` |
| `Cells` | Number of cells (for reference) | `1,108` |
| `Platform` | Sequencing platform | `10x Genomics` |
| `PMID` | PubMed ID of source publication | `31061481` |
| `major.cts` | Major cell types present | `B, CD4Tconv, CD8T, Mono/Macro, NK` |
| `cancer.type` | Cancer type abbreviation | `NSCLC` |

#### Example Metadata Entry

```csv
Dataset.Name,ID,Species,Treatment,Patients,Cells,Platform,PMID,major.cts,cancer.type
NSCLC_GSE127471,T010059,Human,None,1,"1,108",10x Genomics,31061481,"B, CD4Tconv, CD8T, Mono/Macro, NK",NSCLC
Glioma_GSE84465,T010038,Human,None,4,"3,533",Smart-seq2,29091775,"AC-like Malignant, Astrocyte, Mono/Macro, Neuron",Glioma
```

#### Selecting a Dataset to Process

In `pipeline.R`, specify which dataset to process by index:

```r
# Option 1: Hardcoded index
input <- 1  # Process first dataset in metadata CSV

# Option 2: Command line argument (for batch processing)
args <- commandArgs(trailingOnly = TRUE)
input <- as.integer(args[1])

# Get dataset name from metadata
dataset <- meta.tisch[[input, "Dataset.Name"]]
```

### Creating Your Own Metadata File

If you have custom datasets (not from TISCH), create a metadata CSV with at minimum:

```csv
Dataset.Name,Species,Cells
my_dataset_1,Human,5000
my_dataset_2,Mouse,3000
```

Then ensure your data is organized to match:
```
<data.path>/
├── my_dataset_1/
│   ├── *.h5
│   └── *_CellMetainfo_table.tsv
└── my_dataset_2/
    ├── *.h5
    └── *_CellMetainfo_table.tsv
```

### Downloading TISCH Data

To download data from TISCH (Tumor Immune Single-cell Hub):

1. Visit [TISCH2](http://tisch.comp-genomics.org/)
2. Select your dataset of interest
3. Download the expression matrix (H5) and cell metadata (TSV)
4. Place files in the appropriate folder structure

## Quick Start

### Step 1: Prepare Your Data

Ensure your scRNA-seq data is in one of these formats:
- 10x Genomics HDF5 (`.h5`)
- 10x Genomics matrix directory (`filtered_feature_bc_matrix/`)
- TISCH-formatted data with metadata TSV

### Step 2: Configure the Pipeline

Edit `code/pipeline.R` to set your data path and dataset:

```r
# Set your data path
data.path <- "/path/to/your/data"

# Set dataset index or name
input <- 1  # or use command line args
dataset <- "your_dataset_name"
```

### Step 3: Run the Pipeline

From R, in the `sc_rna_pipeline` directory:

```r
# Set working directory
setwd("pipelines/sc_rna_pipeline")

# Run the pipeline
source("code/pipeline.R")
```

Or via command line with SLURM:

```bash
Rscript code/pipeline.R 1  # Pass dataset index as argument
```

## Pipeline Steps

### First Stage: All Cells

| Step | Function | Output | Description |
|------|----------|--------|-------------|
| 1 | `create.seurat.from.Tisch.h5()` | `seurat.raw.RData` | Load counts and metadata |
| 2 | `remove_umap_from_tisch()` | `seurat.raw.clean.RData` | Clean precomputed embeddings |
| 3 | `do.singler()` | `singler.RData` | Cell type predictions |
| 4 | `qc.seurat()` | `seurat.qc.RData` | Filter cells (MT%, features, counts) |
| 5 | `seurat.SCT()` | `seurat.sct.RData` | SCTransform normalization |
| 6 | `seurat.PCA()` | `seurat.pca.RData` | PCA (up to 200 components) |
| 7 | `seurat.findNeighbors()` | `seurat.neighbors.RData` | Build neighbor graph |
| 8 | `seurat.findClusters()` | `seurat.clusters.RData` | Initial clustering |
| 9 | `seurat.visualization()` | `seurat.viz.RData` | UMAP and t-SNE |
| 10 | `do.clusterings()` | `seurat.clusterings.RData` | Multi-resolution clustering |
| 11 | `find.markers.by.all.clusterings()` | `all.markers.RData` | Marker genes per cluster |
| 12 | `get.cell.and.cluster.stats.for.all.clusterings()` | `cell.and.cluster.stats.RData` | Cluster statistics |
| 13 | `do.curated.cell.types()` | `curated.cell.types.RData` | Assign cell type labels |
| 14 | `get.covs()` | `covs.RData` | Combine all covariates |
| 15 | `create.seurat.list.for.ISCVAM()` | `seurat.list.RData` | Package for ISCVAM |

### Second Stage: Per-Layer Analysis

For each major cell type with ≥100 cells:

| Step | Function | Output | Description |
|------|----------|--------|-------------|
| 1 | `second_stage_analyze()` | `<layer>_seurat.list.RData` | Subset and re-analyze |
| 2-10 | (Same as first stage) | Per-layer files | Layer-specific analysis |

### Final Export

| Step | Function | Output | Description |
|------|----------|--------|-------------|
| 1 | `do.multiStages()` | `seurat.list.multi.stages.RData` | Combine all layers |
| 2 | `create.artifacts.for.all.layers()` | `all.layers.RData` | Build heatmap artifacts |
| 3 | `writing.iscvam.h5()` | `<dataset>.test.pipeline.h5` | Export H5 file |

## Output Files

### Intermediate Files (RData)

All intermediate results are saved as `.RData` files in the dataset results folder, allowing resumption from any checkpoint.

### Excel Reports

| File | Description |
|------|-------------|
| `all.markers.xlsx` | Marker genes for all clusterings |
| `cells.and.clusterings.xlsx` | Cell and cluster statistics |
| `<layer>_all.markers.xlsx` | Per-layer marker genes |
| `<layer>_cells.and.clusterings.xlsx` | Per-layer statistics |

### Final H5 File

The output H5 file (`<dataset>.test.pipeline.h5`) contains:

- **assaysData/RNA/**: Normalized expression matrix
- **artifacts/{layer}/covs/**: Cell metadata and embeddings
- **artifacts/{layer}/clusterings/**: Cluster assignments and heatmap data

Structure:
```
<dataset>.test.pipeline.h5
├── assaysData/
│   └── RNA/
│       └── matrix/
├── artifacts/
│   ├── all/
│   │   ├── covs
│   │   ├── discreteCovs
│   │   ├── continuousCovs
│   │   └── clusterings/
│   │       └── clusters_SCT_0.1/
│   │           ├── markers
│   │           ├── all_markers_lognorm
│   │           ├── all_markers_scaled
│   │           └── ...
│   ├── <layer1>/
│   └── <layer2>/
```

## Customization

### QC Parameters

Edit `scRNA_functions_pipeline.R` to adjust QC thresholds:

```r
# In qc.seurat() function
percent.mt <= 20        # Maximum mitochondrial percentage
nFeature_RNA >= 500     # Minimum genes per cell
nCount_RNA >= 1000      # Minimum UMIs per cell
```

### Clustering Resolutions

Default resolutions are adjusted based on dataset size:

- **>5,000 cells**: `c(0.01, 0.05, 0.08, 0.1, 0.2)`
- **≤5,000 cells**: `c(0.01, 0.02, 0.03, 0.05, 0.08, 0.1)`

### PCA Components

Adjust in `pipeline.R`:

```r
seurat.pca <- seurat.PCA(seurat.sct, pcs.to.compute = 200)
seurat.neighbors <- seurat.findNeighbors(seurat.pca, pcs.to.analyze = 40)
```

### Second Stage Threshold

Layers with fewer than 100 cells are skipped. Modify in `scRNA_functions_pipeline.R`:

```r
# Minimum cells for second-stage analysis
min_cells_for_layer <- 100
```

## Troubleshooting

### Memory Issues

For large datasets, increase memory limits:

```r
options(future.globals.maxSize = 65000 * 1024^2)  # ~65GB
```

### Missing SingleR References

Create references if not available:

```r
library(celldex)
refs <- list(
  hpca = HumanPrimaryCellAtlasData(),
  blueprint = BlueprintEncodeData()
)
save(refs, file = "reference_data/singler.refs.rdata")
```

### Dataset Not Found

Verify your data path structure matches expected format:
```
<data.path>/<dataset>/
├── *.h5                        # Count matrix
└── *_CellMetainfo_table.tsv    # Cell metadata (TISCH format)
```

## Support

For detailed pipeline logic, see [pipeline.md](code/pipeline.md).

For issues or questions, please refer to the main ISCVAM documentation or contact the development team.
