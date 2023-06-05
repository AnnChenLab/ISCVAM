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
- Node `v16.20.0` for the backend
- Node `v18.16.0` for the backend

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
