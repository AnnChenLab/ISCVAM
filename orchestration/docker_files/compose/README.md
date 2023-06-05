# Docker Usage for ISCVAM

## Requirements
Be sure to have the following technologies installed with the required version:
 - Docker
	 - With the CLI commands enabled (for running `docker` and `docker-compose`)
	 - https://docs.docker.com/engine/install/

## How to run ISCVAM:
1. Be sure that the port `3000` and `8001` are available to use, also be sure to be in the `orchestration/compose` folder
2. In one terminal run 
```bash
> docker-compose run  
```
- This will download the docker images and mount both server and application

3. You can access to your application in `http://localhost:3000/`
## How to use ISCVAM wit your own data
Before running the docker containers (following the How to run steps) you'll need to:
 1. Generate your ***.h5 file***
 2. Paste your ***.h5 file*** on the *datasets* folder
 3. Add a reference for your file in the ***backend/config/datasets.json*** file:
	 - Let's assume that your file is named ***"test_iscvamR.h5"***
	 - You'll need to add at the end of the ***datasets.json*** the following code:
 ```javascript
{
	// Name: key to identify your dataset
	"name":  "test_iscvamR",
	// Filename: path were your dataset is located
	"filename":  "data/datasets/test_iscvamR.h5"
}
```

 6. Your  ***datasets.json***  will look like this at the end:
  ```javascript
[
	{
		"name":  "discovery_dataset",
		"filename":  "data/demo/test_discover_42.h5"
	},
	{
		"name":  "PBMC_multiome",
		"filename":  "data/demo/PBMC_multiome_2.h5"
	},
	{
		"name":  "internal_validation",
		"filename":  "data/demo/internal_validation_multiome_5.h5"
	},
	// This is what you had to paste in the last step
	{
		"name":  "test_iscvamR",
		"filename":  "data/datasets/test_iscvamR.h5"
	}
]
```
 8. Add a reference for your file in the ***frontend/config/app-settings.json*** file:
	 - You'll need to know the *layers* and *modalities* that are contained in your ***.h5 file***
	 - You'll need to add at the end of the ***app-settings.json*** the following code:
```json
{
	//Any name you want, this is the one that is going to be displayed in the application
	"name":  "TEST VALIDATION",
	// You'll need to change the last part after project/, the key of your dataset, 
	// defined in the previous step
	"url":  "http://localhost:8001/project/test_iscvamR",
	// You'll need to define the layers and modalities that are contained in your .h5 file here
	"layers":  [
		{
			"name":  "RNA",
			"id":  "rna"
		},
		{
			"name":  "ATAC",
			"id":  "atac"
		},
		{
			"name":  "WNN",
			"id":  "wnn"
		}
	],
	"modalities":  [
		"RNA",
		"ATAC"
	]
}
```
 9. Your  ***app-settings.json***  will look like this at the end:
 ```json
{
	"datasets":  [
	{
		"name":  "multiome discovery WNN",
		"url":  "http://localhost:8001/project/discovery_dataset",
		"layers":  [
			{
				"name":  "RNA",
				"id":  "rna"
			},
			{
				"name":  "ATAC",
				"id":  "atac"
			},
			{
				"name":  "WNN",
				"id":  "wnn"
			}
		],
		"modalities":  [
			"RNA",
			"ATAC"
		],
		"zip":  "linkage_discovery.zip"
	},
	{
		"name":  "multiome validation",
		"url":  "http://localhost:8001/project/internal_validation",
		"layers":  [
			{
				"name":  "RNA",
				"id":  "rna"
			},
			{
				"name":  "ATAC",
				"id":  "atac"
			},
			{
				"name":  "WNN",
				"id":  "wnn"
			}
		],
		"modalities":  [
			"RNA",
			"ATAC"
		]
	},
	{
		"name":  "multiome PBMC",
		"url":  "http://localhost:8001/project/PBMC_multiome",
		"layers":  [
			{
				"name":  "RNA",
				"id":  "rna"
			},
			{
				"name":  "ATAC",
				"id":  "atac"
			},
			{
				"name":  "WNN",
				"id":  "wnn"
			}
		],
		"modalities":  [
			"RNA",
			"ATAC"
		]
	},
	{
		"name":  "TEST VALIDATION",
		"url":  "http://localhost:8001/project/test_iscvamR",
		"layers":  [
			{
				"name":  "RNA",
				"id":  "rna"
			},
			{
				"name":  "ATAC",
				"id":  "atac"
			},
			{
				"name":  "WNN",
				"id":  "wnn"
			}
		],
		"modalities":  [
				"RNA",
				"ATAC"
			]
		}
	]
}
``` 
