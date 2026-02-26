First stage (all cells)
Goal: Ingest one TISCH cohort, QC, normalize, reduce dimensions, cluster at multiple resolutions, call markers, annotate, and prepare artifacts for the full dataset.

Dataset selection and setup
Source: pipeline.R
Dataset index (from SLURM args or hardcoded): input <- "1"
meta.tisch <- read.csv("reference_data/meta.data_check.endo.track.progress.man.csv")
Dataset name: dataset <- meta.tisch[[input, "Dataset.Name"]]
Data root: data.path <- ".../TISCH2_download"
Create working dir: res_local_running/<dataset>/, setwd there
H5 output name: <dataset>.test.pipeline.h5
Load data (counts + metadata) from TISCH
create.seurat.from.Tisch.h5(data.path, dataset):
Reads 10x HDF5 counts file: <data.path>/<dataset>/*.h5 via Read10X_h5(...)
Builds a Seurat object with counts (project = dataset)
Reads TSV metadata: <data.path>/<dataset>/*_CellMetainfo_table.tsv
cbinds metadata into seurat@meta.data
Saves: seurat.raw.RData
Clean precomputed UMAP from TISCH metadata
remove_umap_from_tisch(seurat.raw)
Drops columns UMAP_1, UMAP_2
Saves: seurat.raw.clean.RData
Cell-type predictions with SingleR (for later stats/covs)
load("reference_data/singler.refs.rdata") → refs
do.singler(seurat.raw.clean@assays[["RNA"]]@layers[["counts"]], refs)
Produces predictions across multiple panels (“main”/“fine”)
Saves: singler.RData
QC filtering
qc.seurat(seurat.raw.clean, "Human"):
Adds percent.mt (MT genes: ^MT-) and percent.rp (ribosomal: ^RP[LS])
Filters cells: percent.mt <= 20, nFeature_RNA >= 500 (default), nCount_RNA >= 1000
Saves: seurat.qc.RData
Normalization / variance stabilization (SCT)
seurat.SCT(seurat.qc, debatch=FALSE)
Uses SCTransform (glmGamPoi and return.only.var.genes=TRUE for very large datasets)
Saves: seurat.sct.RData
Dimensionality reduction and graph
PCA: seurat.PCA(seurat.sct, pcs.to.compute = 200) with auto-downshift if tiny
Saves: seurat.pca.RData
Neighbors: seurat.findNeighbors(seurat.pca, pcs.to.analyze = 40)
Saves: seurat.neighbors.RData
Clustering and visualization
Clusters: seurat.findClusters(seurat.neighbors)
Saves: seurat.clusters.RData
UMAP (+ TSNE): seurat.visualization(seurat.clusters) with dims = 1:pcs.to.analyze
Saves: seurat.viz.RData
Multi-resolution clusterings (metadata columns)
do.clusterings(seurat.viz, "SCT")
If > 5,000 cells: resolutions c(0.01, 0.05, 0.08, 0.1, 0.2)
Else: c(0.01, 0.02, 0.03, 0.05, 0.08, 0.1)
Adds columns like clusters_SCT_0.1, etc.
Saves: seurat.clusterings.RData
Marker discovery
find.markers.by.all.clusterings(seurat.clusterings)
For each clusters_* column where ≥ 2 groups:
FindAllMarkers with parallelization, computes diff = pct.1 - pct.2
Writes per-clustering CSVs and returns a named list
Saves: all.markers.RData
Excel export: write_all.markers("all.markers.xlsx", all.markers)
Saves: export.rna.markers.RData
Markers for TISCH major cell types (if present)
find.markers.Tisch.major.cts(seurat.clusterings, dataset)
Scans metadata columns matching Celltype.* and runs marker discovery
Writes <dataset>.all.markers.Tisch.major.cts.xlsx
Saves: tisch.markers.RData
Cell and cluster stats (for reporting)
get.cell.and.cluster.stats.for.all.clusterings(seurat.clusterings, singler)
Uses SingleR “blueprint.main” to characterize clusters
Saves: cell.and.cluster.stats.RData
Excel export: write_cells.clusters.stats('cells.and.clusterings.xlsx', ...)
Saves: export.cells.clusters.stats.RData
Curated cell types
do.curated.cell.types(cell.and.cluster.stats)
Assigns cluster-majority labels with cluster IDs
Saves: curated.cell.types.RData
Clinical covariates and combined covs
Clinical/meta covs from TISCH (UMAP columns removed, filtered to QC’d cells):
get.Tisch.covs.clinical(data.path, dataset, seurat.clusterings)
Saves: covs.clinical.RData
Merge covs: get.covs(seurat.clusterings, singler, "Human", covs.clinical, curated.cell.types)
Combines SingleR labels, curated defaults, and clinical/meta
Saves: covs.RData
Pack first-stage results for ISCVAM
create.seurat.list.for.ISCVAM(seurat.clusterings, covs)
Saves: seurat.list.RData
Enrich first-stage with marker sets
update.seurat.list(seurat.list, all.markers, tisch.markers)
Second stage (per-layer re-analysis)
Goal: For each main category with sufficient cells, re-run DR/clustering/markers within that subset to get higher-resolution structure and artifacts per layer; then produce a multi-layer H5.

Determine layers and eligibility
Layers from covs$Celltype..malignancy. → update Seurat Idents
Compute counts per layer; select layers with n ≥ 100
This threshold is enforced via janitor::tabyl + flag do.secondStage
Per-layer analysis
For each eligible layer:
second_stage_analyze(seurat, covs, layer):
Subset Seurat to idents = layer
Remove existing SCT_* metadata columns
Re-run PCA/Neighbors/Clusters/UMAP/TSNE with analyze.seurat_2nd.stage(...)
Adjusts PCs for 100–200 cell ranges (20 PCs/analyzed dims)
do.clusterings(..., "SCT") to add clusters_* resolutions
find.markers.by.all.clusterings(...) and Excel export <layer>_all.markers.xlsx
Stats and Excel export: <layer>_cells.and.clusterings.xlsx
Curated layer cts and get.covs_2ndstage(...) to build layer covs
Return layer bundle: list(seurat = ..., covs = ..., markers = ...)
Saves: <layer>_seurat.list.RData
Combine layers
do.multiStages(seurat.list):
Runs the per-layer loop above
Returns: seurat.list.all.layers <- c(seurat.list, all.second.stages)
Build artifacts for all layers
create.artifacts.for.all.layers(seurat.list.all.layers):
For each layer:
Clustering heatmap artifacts via create.artifacts.ALL.clusterings(...):
Uses heatmap_artifacts_from_seurat_V5(...) to produce:
markers, all_markers_lognorm, all_markers_scaled, heatmap_markers_lognorm, heatmap_markers_scaled
Covariate artifacts via:
prep_covs_Tisch(covs):
Removes QC/clustering columns to extract discrete clinical covs
layer_artifacts_from_seurat(...):
Assembles covs with required columns: tsne_1, tsne_2, umap_1, umap_2, QC features, clusterings (from both Seurat metadata and TISCH covs)
Attach clusterings artifacts into the layer’s covs bundle
Returns: layers (named by layer: all, plus each selected main category)
Write multi-layer H5
writing.iscvam.h5(seurat.list.all.layers, layers, dataset, fn.h5)
Internally: write_h5(fn, seurat, layers, assays=c("RNA"))
Saves counts under assaysData/RNA/matrix/*
Stores per-layer artifacts/{layer}/covs, discreteCovs, continuousCovs
Stores per-layer clustering heatmap artifacts under artifacts/{layer}/clusterings/...
Quick artifact map
First stage main outputs (in res_local_running/<dataset>/)

Seurat checkpoints: seurat.raw*.RData, seurat.qc.RData, seurat.sct.RData, seurat.pca.RData, seurat.neighbors.RData, seurat.clusters.RData, seurat.viz.RData, seurat.clusterings.RData
SingleR: singler.RData
Markers: all.markers.RData, all.markers.xlsx, <dataset>.all.markers.Tisch.major.cts.xlsx
Stats: cell.and.cluster.stats.RData, cells.and.clusterings.xlsx, export.cells.clusters.stats.RData
Covariates/annotations: covs.clinical.RData, curated.cell.types.RData, covs.RData
ISCVAM bundle: seurat.list.RData
Second stage main outputs

Per-layer bundles saved as <layer>_seurat.list.RData
Layer marker and stats Excel files: <layer>_all.markers.xlsx, <layer>_cells.and.clusterings.xlsx
Combined layers: seurat.list.multi.stages.RData, all.layers.RData
Final H5: <dataset>.test.pipeline.h5
If you’d like, I can add a one-page README in res_local_running/<dataset>/ summarizing these steps and artifacts for quick reference, or generate a simple driver that toggles first vs. second stage for batch runs.