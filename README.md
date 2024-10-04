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

## Add new dataset

1. **Upload the dataset** at: "/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/data"
2. Update 3 following file:

   1. Dataset name and location in the "datasets.json"
      go to "/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/backend/datasets.json"
      insert the **dataset name** and the **link** to the data
      **Example:¬**

   ```bash
     {
         "name": "ABC_dataset",
         "filename": "/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/data/ABC_dataset.h5"
     }
   ```
   <br/>

   2. Frontend setting in the "frontend/build/app-settings.json"
      go to "/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/frontend/build/app-settings.json"
      insert the following at the end of the file inside the datasets: []
      **Example:**

   ```bash
     {
           "name": "ABC-example", 
           "url": "https://chenlab.chpc.utah.edu/backend/project/ABC-example", 
           "layers": [
               {
                   "name": "all", 
                   "id": "all"
               }

           ], 
           "modalities":[
               "RNA"
           ]
       }
   ```
   <br/>

   3. Public frontend setting in the "frontend/public/app-settings.json"
      go to "/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/frontend/public/app-settings.json"
      insert the same config from 2.2 at the end of the file inside the datasets: []
      **Example:**

   ```bash
     {
           "name": "ABC-example", 
           "url": "https://chenlab.chpc.utah.edu/backend/project/ABC-example", 
           "layers": [
               {
                   "name": "all", 
                   "id": "all"
               }

           ], 
           "modalities":[
               "RNA"
           ]
     }
   ```
