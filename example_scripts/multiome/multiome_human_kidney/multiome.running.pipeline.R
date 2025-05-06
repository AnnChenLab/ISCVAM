

data.path <- "/project/multiome/multiome_human_kidney/raw_data"

source("/project/multiome/pipeline.standard.codes/updated_write_h5_scripts_using_iscvamR.R")
source("/project/multiome/pipeline.standard.codes/multiome_functions_pipeline.R")
source("/project/multiome/pipeline.standard.codes/iscvamR_h5_customized.SeuratV5.R")

filename <- file.path(data.path, "multiome_human_kidney_update.res.h5")
path.result <- "/project/multiome/multiome_human_kidney/results_update.resolutions"

load("/project/multiome/pipeline.standard.codes/singler.refs.rdata")
species <- "Human"

options(future.globals.maxSize = 50 * 1024^3) #here 50GB

##########
#creating seurat object
seurat.raw <- create_seurat_multiome(data.path, meta.data=FALSE)
save(seurat.raw, file = file.path(path.result, "seurat.raw.RData"))

singler <- do.singler(seurat.raw, refs)
save(singler, file = file.path(path.result, "singler.RData"))

peaks.calling <- call.peaks_from_signac(seurat.raw, data.path)
save(peaks.calling, file = "peaks.calling.RData")

seurat.peaks <- create_peaks_assay(peaks.calling, seurat.raw)
save(seurat.peaks, file = "seurat.peaks.RData")




dimension.before.qc <- retrieve_dimension(seurat.peaks)
print(dimension.before.qc)

############
#QC steps
seurat.qc_rna <- qc_RNA(seurat.peaks, "Human")
save(seurat.qc_rna, file = file.path(path.result,"seurat.qc_rna.RData"))

seurat.qc_atac <- qc_ATAC(seurat.qc_rna, meta.data =FALSE)
save(seurat.qc_atac, file = file.path(path.result,"seurat.qc_atac.RData"))


dimension.after.qc <- retrieve_dimension(seurat.qc_rna)

############
#clustering steps for both RNA and ATAC

seurat.analyzed.rna <- analyze_rna(seurat.qc_atac, debatch=FALSE)
save(seurat.analyzed.rna, file = file.path(path.result,"seurat.analyzed.rna.RData"))

clustering.rna <- do.clusterings(seurat.analyzed.rna, "SCT", resolution.lst)
save(clustering.rna, file = file.path(path.result,"clustering.rna.RData"))

seurat.analyzed.atac <- analyze_atac(clustering.rna, debatch=FALSE)
save(seurat.analyzed.atac, file = file.path(path.result,"seurat.analyzed.atac.RData"))

clustering.atac <- do.clusterings(seurat.analyzed.atac, "peaks", resolution.lst)
save(clustering.atac, file = file.path(path.result,"clustering.atac.RData"))

seurat.wnn <- build_WNN(clustering.atac, debatch.atac=FALSE)
save(seurat.wnn, file = file.path(path.result,"seurat.wnn.RData"))

seurat.wnn.clusterings <- do.clusterings.WNN(seurat.wnn, resolution.lst)
save(seurat.wnn.clusterings, file = file.path(path.result,"seurat.wnn.clusterings.RData"))

##############
#generating markers for all clusterings, each clustering will have 2 markers for each assay ("SCT" and "peaks")

markers_sct_clusterings <- find.markers.by.all.clusterings.for.two.assays(seurat.wnn.clusterings, group=0, "SCT", resolution.lst)                           
markers_sct_clusterings <- formating_markers_clusterings(markers_sct_clusterings)
save(markers_sct_clusterings, file = file.path(path.result,"markers_sct_clusterings.RData"))

markers_peaks_clusterings <- find.markers.by.all.clusterings.for.two.assays(seurat.wnn.clusterings, group=0, "peaks", resolution.lst)
markers_peaks_clusterings <- formating_markers_clusterings(markers_peaks_clusterings)
save(markers_peaks_clusterings, file = file.path(path.result,"markers_peaks_clusterings.RData"))

markers_wnn_clusterings <- find.markers.by.all.clusterings.for.two.assays(seurat.wnn.clusterings, group=1, "wnn", resolution.lst)   
markers_wnn_clusterings <- formating_markers_clusterings(markers_wnn_clusterings)
save(markers_wnn_clusterings, file = file.path(path.result,"markers_wnn_clusterings.RData"))


##########################

cell.and.cluster.stats <- get.cell.and.cluster.stats.for.all.clusterings(seurat.wnn.clusterings, singler)
save(cell.and.cluster.stats, file = file.path(path.result,"cell.and.cluster.stats.RData"))

export.cells.clusters.stats <- write_cells.clusters.stats('cells.and.clusterings.xlsx', cell.and.cluster.stats, seurat.wnn.clusterings)

default.cts <- do.curated.cell.types(cell.and.cluster.stats)
save(default.cts, file = file.path(path.result,"default.cts.RData"))

seurat.wnn.clusterings$default.cts <- default.cts

covs <- get.covs(seurat.wnn.clusterings, singler, species, default.cts) 
save(covs, file = file.path(path.result,"covs.RData"))

###############

markers_sct_clusterings_update <- default.cts_find.markers.and.add.to.markers.sct.clusterings(seurat.wnn.clusterings,
                                                                                              default.cts, 
                                                                                              markers_sct_clusterings)

markers_sct_clusterings_update <- formating_markers_clusterings(markers_sct_clusterings_update)
save(markers_sct_clusterings_update, file = file.path(path.result,"markers_sct_clusterings_update.RData"))

###############
#prep for ISCVAM-h5 compliant

seurat.wnn.clusterings$default.cts <- default.cts
save(seurat.wnn.clusterings, file = file.path(path.result,"seurat.wnn.clusterings.RData"))


clus.avg.for.sct.clusterings <- create.artifacts.ALL.clusterings(seurat.wnn.clusterings, markers_sct_clusterings_update)
save(clus.avg.for.sct.clusterings, file = file.path(path.result,"clus.avg.for.sct.clusterings.RData"))

clus.avg.for.peaks.clusterings <- create.artifacts.ALL.clusterings(seurat.wnn.clusterings, markers_peaks_clusterings)
save(clus.avg.for.peaks.clusterings, file = file.path(path.result,"clus.avg.for.peaks.clusterings.RData"))

clus.avg.for.wnn.clusterings <- create.artifacts.ALL.clusterings(seurat.wnn.clusterings, markers_wnn_clusterings)
save(clus.avg.for.wnn.clusterings, file = file.path(path.result,"clus.avg.for.wnn.clusterings.RData"))

################
heatmap_artifacts <- assemble_heatmap_artifacts(clus.avg.for.sct.clusterings, 
                                                clus.avg.for.peaks.clusterings, 
                                                clus.avg.for.wnn.clusterings)

save(heatmap_artifacts, file = file.path(path.result,"heatmap_artifacts.RData"))


############
#add gene names to counts matrix
path.10x <- paste0(data.path, "/", list.files(data.path, pattern = ".*.h5"))
counts = Read10X_h5(path.10x)
seurat.wnn.clusterings@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]] <- counts[["Gene Expression"]]@Dimnames[[1]]
seurat.wnn.clusterings@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]] <- counts[["Gene Expression"]]@Dimnames[[2]]

save(seurat.wnn.clusterings, file = file.path(path.result,"seurat.wnn.clusterings.RData"))


#############
#writing iscvam h5-compliant
write_mm_h5(seurat.wnn.clusterings, covs, heatmap_artifacts, filename) 
##########



