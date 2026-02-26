# ISCVAM - Interactive Single Cell Visual Analytics for Multiomics

## About

**ISCVAM** is a fast, interactive tool for visualizing and investigating *single-cell multi-omics* data. This repository contains:

1. **Web Application** - Frontend and backend code for the ISCVAM visualization platform
2. **Analysis Pipelines** - R pipelines for processing single-cell data into ISCVAM-compatible H5 format

ISCVAM can be accessed at https://chenlab.utah.edu/iscvam/

Contact: ann.chen@hci.utah.edu

---

## Analysis Pipelines

We provide two analysis pipelines for preparing data for ISCVAM:

### 1. scRNA-seq Pipeline
For single-cell RNA sequencing data only.
- Location: `pipelines/sc_rna_pipeline/`
- Input: 10x Genomics scRNA-seq data
- Output: H5 file with expression data, clustering, and annotations
- See [scRNA-seq Pipeline README](pipelines/sc_rna_pipeline/README.md) for detailed instructions

### 2. Multiome Pipeline
For 10x Multiome data (scRNA-seq + scATAC-seq).
- Location: `pipelines/multiome_pipeline/`
- Input: 10x Genomics Multiome data (Cell Ranger ARC output)
- Output: H5 file with RNA, ATAC, and integrated WNN analysis
- See [Multiome Pipeline README](pipelines/multiome_pipeline/README.md) for detailed instructions

---

## Web Application

### Requirements

Be sure to have the following technologies installed with the required version:

- Docker
  - With the CLI commands enabled (for running `docker` and `docker-compose`)
  - https://docs.docker.com/engine/install/
- If not using Docker, make sure you have:
  - For the backend:
    - Node `v16.20.0`
    - node hdf5 addon: https://github.com/zhihua-chen/hdf5.node
  - For the frontend:
    - Node `v18.16.0`

### Folder Structure

```bash
ISCVAM/
├── backend/                    # Backend server (Node.js)
├── frontend/                   # Frontend web application (React)
├── pipelines/
│   ├── sc_rna_pipeline/        # scRNA-seq analysis pipeline
│   └── multiome_pipeline/      # Multiome (RNA + ATAC) analysis pipeline
├── example_scripts/            # Example processing scripts
│   ├── multiome/
│   └── scRNA_seq/
├── orchestration/
│   └── docker_files/
│       ├── backend/
│       ├── frontend/
│       ├── pipeline/
│       └── compose/            # Docker Compose configuration
│           ├── backend/config/datasets.json
│           ├── frontend/config/app-settings.json
│           ├── datasets/       # Place your .h5 files here
│           └── docker-compose.yml
└── 
```

---

## Example Datasets

We applied ISCVAM to investigate cell populations using multiple *multiome* datasets for proof of principle:

| Dataset | Cells | Description |
|---------|-------|-------------|
| Ovarian Cancer TRM/ReCir | 16,232 | CD8+ T cells from 4 patients (Anadon, Yu et al. 2022, GEO: GSE192780) |
| PBMC | 11,172 | 10x Genomics PBMC sample |
| Human Brain | 3,233 | 10x Genomics healthy brain tissue |
| Human Kidney Cancer | 22,722 | 10x Genomics kidney cancer nuclei |

In our *manuscript*, datasets were organized as:
- **Discovery dataset**: Paired TRM and ReCir samples from patient 100809M (ovarian cancer)
- **Internal validation dataset**: Remaining 3 patients with paired TRM/ReCir samples
- **External validation dataset**: PBMC multiome data

### References
- Anadon et al. *Ovarian cancer immunogenicity is governed by a narrow subset of progenitor tissue-resident memory T cell.* Cancer Cell (2022)
- [10x Multiome PBMC](https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-1-0-0)
- [10x Multiome Human Brain](https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0)
- [10x Multiome Human Kidney Cancer](https://www.10xgenomics.com/datasets/human-kidney-cancer-nuclei-isolated-with-chromium-nuclei-isolation-kit-saltyez-protocol-and-10x-complex-tissue-dp-ct-sorted-and-ct-unsorted-1-standard) 

