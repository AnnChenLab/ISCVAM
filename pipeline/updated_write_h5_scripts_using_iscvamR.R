calculate_gene_activity <- function(seurat) {
  # # Set annotation for peaks in the ATAC assay
  # Annotation(seurat@assays$peaks) <- seurat@assays$ATAC@annotation
  # 
  # library(EnsDb.Hsapiens.v86)
  # # get gene annotations for hg38
  # annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  # seqlevelsStyle(annotation) <- "UCSC"
  
  
  # Compute gene activities
  gene_activities <- GeneActivity(seurat, assay = "peaks")
  
  # Add gene activities to Seurat object
  seurat[["GeneActivity"]] <- CreateAssayObject(counts = gene_activities)
  
  return(seurat)
}


attach_closest_features_to_atac_markers <- function(seurat, heatmap_artifacts) {
  # # Setting annotation for ATAC assay
  # Annotation(seurat@assays$peaks) <- seurat@assays$ATAC@annotation
  
  # Calculate closest feature
  closest_feature <- ClosestFeature(seurat@assays$peaks, regions = granges(seurat@assays$peaks))
  
  # Annotate markers function
  annotate_markers <- function(markers) {
    markers %>%
      left_join(closest_feature[, c("query_region", "gene_name", "type", "distance")],
                by = c("gene" = "query_region")) %>%
      select(gene, gene_name, type, distance, everything()) %>%
      mutate_if(is.factor, as.character)
  }
  
  # Iterate over each layer and clustering to annotate markers
  for (layer in names(heatmap_artifacts)) {
    print(paste0('working on ', layer))
    for (clustering in names(heatmap_artifacts[[layer]]$clusterings)) {
      current_clustering <- heatmap_artifacts[[layer]]$clusterings[[clustering]]
      print(paste0('   clustering: ', clustering))
      # If markers exist in the current clustering and layer is ATAC, annotate markers
      if ("markers" %in% names(current_clustering) && layer == "atac") {
        print('annotating single modal atac markers')
        heatmap_artifacts[[layer]]$clusterings[[clustering]]$markers <- annotate_markers(current_clustering$markers)
      } 
      
      # If ATAC is present in the current clustering, annotate markers
      else if ("atac" %in% names(current_clustering)) {
        heatmap_artifacts[[layer]]$clusterings[[clustering]]$atac$markers <- annotate_markers(current_clustering$atac$markers)
      }
    }
  }
  return(heatmap_artifacts)
}


# assemble_heatmap_artifacts <- function(clus.avg.ALL.layers_h5,
#                                        clus.avg_main.cat_h5,
#                                        clus.avg_sub.cat_h5,
#                                        clus.avg_default.ct_h5,
#                                        clus.avg_wnn_0.6_h5,
#                                        clus.avg_curated.ct_h5) {
#   list(
#     rna = c(
#       clus.avg.ALL.layers_h5$rna,
#       list(
#         main.categories=clus.avg_main.cat_h5,
#         sub.categories=clus.avg_sub.cat_h5,
#         default.ct=clus.avg_default.ct_h5
#       )),
#     atac = clus.avg.ALL.layers_h5$atac,
#     wnn = list(
#       curated.ct_wnn=clus.avg_curated.ct_h5,
#       clusters_wsnn_0.6=clus.avg_wnn_0.6_h5
#     )
#   )
# }

# heatmap_artifacts = assemble_heatmap_artifacts(clus.avg.ALL.layers_h5,
#                                                clus.avg_main.cat_h5,
#                                                clus.avg_sub.cat_h5,
#                                                clus.avg_default.ct_h5,
#                                                clus.avg_wnn_0.6_h5,
#                                                clus.avg_curated.ct_h5)

write_mm_h5 <- function(seurat, covs_discrete, heatmap_artifacts, filename) {
  
  # Calculate gene activity and attach closest features
  seurat <- seurat %>% calculate_gene_activity()
  heatmap_artifacts <- attach_closest_features_to_atac_markers(seurat, heatmap_artifacts)
  
  # Define QC features and clustering features
  qc_features <- c("nFeature_RNA", "nCount_RNA", "nFeature_ATAC", "nCount_ATAC",
                   "nFeature_peaks", "nCount_peaks", "percent.mt", "percent.rp",
                   "nucleosome_signal", "nucleosome_percentile", "TSS.enrichment", "TSS.percentile")
  
  clustering_features <- grep("clusters.*", colnames(seurat@meta.data), value = TRUE)
  
  # Define reductions
  reductions <- list(rna = c("tsne.rna", "umap.rna"),
                     atac = c("tsne.atac", "umap.atac"),
                     wnn = c("wnn.tsne", "wnn.umap"))
  
  # Process each layer and attach heatmap artifacts
  layer_names <- c("rna", "atac", "wnn")
  layers <- lapply(layer_names, function(layer) {
    layer_artifacts <- layer_artifacts_from_seurat(
      seurat, 
      qc_features = qc_features, 
      umap = reductions[[layer]][2],
      tsne = reductions[[layer]][1],
      clustering_names = clustering_features,
      extra_discrete_covs = covs_discrete
    )
    layer_artifacts[["clusterings"]] <- heatmap_artifacts[[layer]]
    layer_artifacts
  })
  
  names(layers) <- layer_names
  
  # Write to .h5 file
  write_h5(filename, seurat, layers = layers, assays = c("RNA", "peaks", "GeneActivity"))
}
