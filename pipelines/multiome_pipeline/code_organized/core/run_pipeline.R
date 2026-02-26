# ============================================================================
# MULTIOME PIPELINE RUNNER
# ============================================================================
# This script runs the complete multiome analysis pipeline
# Make sure to run configure_dataset.R first to set your dataset

# Check if configuration is loaded
if (!exists("dataset_name") || !exists("species")) {
  cat("❌ ERROR: Dataset configuration not found!\n")
  cat("Please run: source('code_organized/configure_dataset.R') first\n")
  stop("Configuration required")
}

cat("=== STARTING MULTIOME PIPELINE ===\n")
cat("Dataset:", dataset_name, "\n")
cat("Species:", species, "\n")
cat("===============================\n\n")

# Load required functions
cat("Loading required functions...\n")
#load("code_organized/ref_data/singler.refs.rdata")

source("code_organized/utils/iscvamR_h5_customized.SeuratV5.R")
source("code_organized/utils/updated_write_h5_scripts_using_iscvamR.R")
source("code_organized/core/multiome_functions_pipeline.R")

# Create seurat object
cat("Step 1: Creating Seurat object from:", data.path, "\n")
seurat.raw <- create_seurat_multiome(data.path, meta.data=FALSE)

# Change to results directory
cat("Step 2: Setting up results directory...\n")
data.output <- results_dir
save(seurat.raw, file = file.path(data.output, "seurat.raw.RData"))

cat("Step 3: Running SingleR cell type annotation...\n")

singler <- do.singler(seurat.raw, refs)
save(singler, file = file.path(data.output,"singler.RData"))

cat("Step 4: Peak calling with MACS3...\n")
peaks.calling <- call.peaks_from_signac(seurat.raw, data.path)
seurat.peaks <- create_peaks_assay(peaks.calling, seurat.raw)
save(seurat.peaks, file = file.path(data.output,"seurat.peaks.RData"))

dimension.before.qc <- retrieve_dimension(seurat.raw)

cat("Step 5: Quality control...\n")

seurat.qc_rna <- qc_RNA(seurat.peaks, species)
save(seurat.qc_rna, file = file.path(data.output,"seurat.qc_rna.RData"))

seurat.qc_atac <- qc_ATAC(seurat.qc_rna, meta.data=FALSE)
save(seurat.qc_atac, file = file.path(data.output,"seurat.qc_atac.RData"))

dimension.after.qc <- retrieve_dimension(seurat.qc_rna)


#cat("Step 5.1: Gene Activity... \n") # (Optional)

# seurat.gene.activities <- calculate_gene_activity(seurat.qc_rna)
# save(seurat.gene.activities, file = file.path(data.output,"seurat.gene.activities.RData"))


cat("Step 6: RNA analysis and clustering...\n")
seurat.analyzed.rna <- analyze_rna(seurat.qc_atac, debatch=FALSE)
save(seurat.analyzed.rna, file = file.path(data.output,"seurat.analyzed.rna.RData"))

clustering.rna <- do.clusterings(seurat.analyzed.rna, "SCT", resolution.lst)
save(clustering.rna, file = file.path(data.output,"clustering.rna.RData"))

cat("Step 7: ATAC analysis and clustering...\n")
seurat.analyzed.atac <- analyze_atac(clustering.rna, debatch=FALSE)
save(seurat.analyzed.atac, file = file.path(data.output,"seurat.analyzed.atac.RData"))

clustering.atac <- do.clusterings(seurat.analyzed.atac, "peaks", resolution.lst)
save(clustering.atac, file = file.path(data.output,"clustering.atac.RData"))

cat("Step 8: WNN integration...\n")
seurat.wnn <- build_WNN(clustering.atac, debatch.atac=FALSE)
save(seurat.wnn, file = file.path(data.output,"seurat.wnn.RData"))

seurat.wnn.clusterings <- do.clusterings.WNN(seurat.wnn, resolution.lst)
save(seurat.wnn.clusterings, file = file.path(data.output,"seurat.wnn.clusterings.RData"))

# ============================================================================
# Step 8.1: Azimuth Annotation
# ============================================================================
cat("\nStep 8.1: Azimuth Annotation...\n")
# one-time (downloads the PBMC Azimuth reference)
#check available dataset in seurat run
#AvailableData()[grep("kidney|motor|cortex|ref", AvailableData()$Dataset, ignore.case=TRUE), c("Dataset","Version")]
# to install choose a dataset run
#SeuratData::InstallData("pbmcref")
# SeuratData::InstallData("kidneyref")
# SeuratData::InstallData("humancortexref")
# seurat.wnn.clusterings <- azimuth_annotation(
#                               seurat.wnn.clusterings,
#                               reference = "humancortexref",
#                               pred_cols = c(
#                                 "predicted.class",
#                                 "predicted.subclass",
#                                 "predicted.cluster",
#                                 "predicted.cross_species_cluster"
#                               ),
#                               join_layers = FALSE
#                             )
# seurat.wnn.clusterings <- azimuth_annotation(
#                             seurat.wnn.clusterings,
#                             reference = "kidneyref",
#                             pred_cols = c("predicted.celltype.l1", "predicted.celltype.l2", "predicted.celltype.l3"),
#                             join_layers = FALSE
#                           )
seurat.wnn.clusterings <- azimuth_annotation(
  seurat.wnn.clusterings,
  reference = "kidneyref",
  pred_cols = c("predicted.annotation.l1", "predicted.annotation.l2", "predicted.annotation.l3"),
  join_layers = FALSE
)


save(seurat.wnn.clusterings, file = file.path(results_dir, "seurat.wnn.clusterings.RData"))

# p1<-DimPlot(
#   seurat.wnn.clusterings,
#   reduction = "umap.rna",
#   group.by = "azimuth_kidneyref_predicted.annotation.l1",
#   label = TRUE
# )
# p2<-DimPlot(
#   seurat.wnn.clusterings,
#   reduction = "umap.rna",
#   group.by = "azimuth_kidneyref_predicted.annotation.l2",
#   label = TRUE
# )
# p1|p2

cat("Step 9: Finding markers...\n")
markers_sct_clusterings <- find.markers.by.all.clusterings.for.two.assays(seurat.wnn.clusterings, group=0, "SCT", resolution.lst)                           
markers_sct_clusterings <- formating_markers_clusterings(markers_sct_clusterings)
save(markers_sct_clusterings, file = file.path(data.output,"markers_sct_clusterings.RData"))

markers_peaks_clusterings <- find.markers.by.all.clusterings.for.two.assays(seurat.wnn.clusterings, group=0, "peaks",resolution.lst )
markers_peaks_clusterings <- formating_markers_clusterings(markers_peaks_clusterings)
save(markers_peaks_clusterings, file = file.path(data.output,"markers_peaks_clusterings.RData"))

markers_wnn_clusterings <- find.markers.by.all.clusterings.for.two.assays(seurat.wnn.clusterings, group=1, "wnn", resolution.lst)   
markers_wnn_clusterings <- formating_markers_clusterings(markers_wnn_clusterings)
save(markers_wnn_clusterings, file = file.path(data.output,"markers_wnn_clusterings.RData"))

cat("Step 10: Cell type annotation and statistics...\n")
cell.and.cluster.stats <- get.cell.and.cluster.stats.for.all.clusterings(seurat.wnn.clusterings, singler)
save(cell.and.cluster.stats, file = file.path(data.output,"cell.and.cluster.stats.RData"))

export.cells.clusters.stats <- write_cells.clusters.stats('cells.and.clusterings.xlsx', cell.and.cluster.stats, seurat.wnn.clusterings)

# default.cts <- do.curated.cell.types2(
#   cell.and.cluster.stats,
#   mode = "majority",
#   clustering = "clusters_wsnn_0.8"
# )
default.cts <- do.curated.cell.types2(
                                    cell.and.cluster.stats,
                                    mode = "direct",
                                    direct_label_col = "azimuth_kidneyref_predicted.annotation.l1",
                                    append_cluster = FALSE
                                  )

save(default.cts, file = file.path(data.output,"default.cts.RData"))

# Generate covs with all metadata
covs <- get.covs(seurat.wnn.clusterings, singler, species, default.cts)

cat("  Covs columns:", paste(colnames(covs), collapse=", "), "\n") 
save(covs, file = file.path(data.output,"covs.RData"))



###############
cat("Step 11: Assigning a cluster resolution (from manual review) as the default cluster...\n")

#adding markers with 2 assays for default.cts
## 
out <- default.cts_find.markers.and.add(
                                seurat.wnn.clusterings = seurat.wnn.clusterings,
                                default.cts = default.cts,
                                markers_lists = list(
                                  sct = markers_sct_clusterings,
                                  peaks = markers_peaks_clusterings
                                ),
                                resolution.lst = resolution.lst,
                                save_file = "markers_default.cts_2assays.RData"
                              )
##
markers_sct_clusterings   <- out$sct
markers_peaks_clusterings <- out$peaks
##
markers_sct_clusterings <- formating_markers_clusterings(markers_sct_clusterings)
markers_peaks_clusterings <- formating_markers_clusterings(markers_peaks_clusterings)

#reorder clusters
order_default.cts <- annotate.curated.cell.types_2assays(default.cts, markers_sct_clusterings)

# Reorder clusters based on a desired sequence
seurat.wnn.clusterings$default.cts <- default.cts

seurat.wnn.clusterings$default.cts <- factor(seurat.wnn.clusterings$default.cts,
                                             levels= order_default.cts$cts.order)
save(seurat.wnn.clusterings, file = file.path(data.output,"seurat.wnn.clusterings.RData"))

markers_sct_clusterings$default.cts <- order_default.cts$markers
save(markers_sct_clusterings, file = file.path(data.output,"markers_sct_clusterings.RData"))
markers_peaks_clusterings$default.cts <- order_default.cts$markers
save(markers_peaks_clusterings, file = file.path(data.output,"markers_peaks_clusterings.RData"))
markers_wnn_clusterings$default.cts <- order_default.cts$markers
save(markers_wnn_clusterings, file = file.path(data.output,"markers_wnn_clusterings.RData"))


cat("Step 12: Preparing ISCVAM artifacts...\n")
clus.avg.for.sct.clusterings <- create.artifacts.ALL.clusterings(seurat.wnn.clusterings, markers_sct_clusterings)
save(clus.avg.for.sct.clusterings, file = file.path(data.output,"clus.avg.for.sct.clusterings.RData"))

clus.avg.for.peaks.clusterings <- create.artifacts.ALL.clusterings(seurat.wnn.clusterings, markers_peaks_clusterings)
save(clus.avg.for.peaks.clusterings, file = file.path(data.output,"clus.avg.for.peaks.clusterings.RData"))

clus.avg.for.wnn.clusterings <- create.artifacts.ALL.clusterings(seurat.wnn.clusterings, markers_wnn_clusterings)
save(clus.avg.for.wnn.clusterings, file = file.path(data.output,"clus.avg.for.wnn.clusterings.RData"))

heatmap_artifacts <- assemble_heatmap_artifacts(clus.avg.for.sct.clusterings, 
                                               clus.avg.for.peaks.clusterings, 
                                               clus.avg.for.wnn.clusterings)
save(heatmap_artifacts, file = file.path(data.output,"heatmap_artifacts.RData"))

cat("Step 13: Exporting final H5 file...\n")
write_mm_h5(seurat.wnn.clusterings, covs, heatmap_artifacts, filename,join_rna=FALSE) 

# Final summary
cat("\n🎉 PIPELINE COMPLETED SUCCESSFULLY! 🎉\n")
cat("========================================\n")
cat("Dataset:", dataset_name, "\n")
cat("Results directory:", results_dir, "\n")
cat("Final H5 file:", filename, "\n")
cat("Full path:", normalizePath(file.path(results_dir, filename)), "\n")
cat("========================================\n")

# Return to main directory
setwd("..")