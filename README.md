# ISCVAM

## About

Scripts for ISCVAM, a tool for exploration of multimodal single cell data.

ISCVAM can be accessed at http://iscva.moffitt.org/

Contact: ann.chen@moffitt.org

# Development

## Requirements

Be sure to have the following technologies installed with the required version:

- Docker
  - With the CLI commands enabled (for running `docker` and `docker-compose`)
  - https://docs.docker.com/engine/install/
- If not using docker, make sure you have
  - for the backend
    - Node `v16.20.0`
    - node hdf5 addon: https://github.com/zhihua-chen/hdf5.node
  - for the frontend
    - Node `v18.16.0`

## Folder Structure:

```bash
    |- backend # backend files
    |- frontend # frontend files
    |- example_scripts
        |- multiome_data
            |- discovery_dataset
            |- internal_validation_dataset
            |- external_validation_dataset
        |- scRNA_seq
            |- TISCH_iscvamR_h5
    |- orchestration
        |- docker-files
            |- backend
            |- frontend
            |- pipeline
            |- compose # docker-compose related files
                |- backend 
                    |- config
                        |- datasets.json
                |- frontend
                    |- config
                        |- app-settings.json
                |- datasets # place to put own .h5 files
                |- docker-compose.yml


    |- pipeline # r files
```

# ISCVAM example datasets

We applied ISCVAM to investigate cell populations using two Multiome datasets for proof of principle. The two multiome datasets are: 
1) 16,232 cells from the matched pairs of sorted CD8+ Tissue resident memory (TRM) and recirculating (ReCir) T cells from 4 ovarian cancer patients in Anadon, Yu et al. 2022
2) 11,172 cells from a PBMC sample from 10X Genomics

In our manucript, the two datasets was used as:
- discovery dataset: which is the paired TRM and ReCir samples from a single patient (patient 100809M) in ovarian cancer paper
- internal validation dataset: the rest of 3 patients with their paired TRM and ReCir samples for each patient from the ovarian cancer paper
- external validation dataset: PBMC multiome data

### Reference datasets
- Anadon et al. Ovarian cancer immunogenicity is governed by a narrow subset of progenitor tissue-resident memory T cell. Cancer Cell (2022) 
- PBMC sample: https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-1-0-0 

