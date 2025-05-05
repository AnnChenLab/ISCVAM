# THIS SCRIPT CUSTOMIZED FOR SEURAT V5 (meaning add @layer[['count']])

#' Generate ISCVAM compliant h5 files.
#' 
#' @description
#' Given a seurat object and various associated artifacts, generate an
#' ISCVAM compliant h5 file.
#' @param fn This is the name of the h5 file.
#' @param seurat This is the seurat object containing the underlying assays.
#' @param layers This is a list of analytical layers for this dataset. Each
#'     layer should contain a matrix called covs, a vector of strings called discreteCovs,
#'     and a vector of strings called continuousCovs. To make use of the heatmap visualizations
#'     in ISCVAM, a list of clustering artifacts should be stored at "clusterings" in this list.
#' @param assays A vector of the name of the assays in the seurat object to be saved to the h5 file.
#' @examples
#' \dontrun{
#'    write_h5("my.iscvam.h5", seurat, layers)
#' }
#' @export
write_h5 <- function(fn, seurat, layers, assays = c("RNA", "peaks", "GeneActivity")) {
  stopifnot(!file.exists(fn))

  rhdf5::h5createFile(fn)
  rhdf5::h5createGroup(fn, "assaysData")

  for (assay in assays) {
    if (assay %in% names(seurat@assays)) {
     if (assay == "RNA"){
       #SeuratV5
       counts.matrix <- seurat@assays[[assay]]@layers[['counts']]
       colnames(counts.matrix) <- colnames(seurat)
      write_h5_extra_assay(fn, counts.matrix, assay)
     } 
      if ((assay == "peaks") | (assay == "GeneActivity")){
       write_h5_extra_assay(fn, seurat@assays[[assay]]@counts, assay)
     } 
    }
  }

  rhdf5::h5createGroup(fn,"artifacts")

  for (artifact_name in names(layers)) {
    #layers <- c("all", "lymphoids", "myloids", "others")
    #or in multiome: layer <- c("rna", "atac", "wnn")
    #to create this list of layers, use function "layer_artifacts_from_seurat" below (output = covs)
    write_h5_artifact(fn, artifact_name, layers[[artifact_name]])

    if (("clusterings" %in% names(layers[[artifact_name]])) & (length(layers[[artifact_name]][["clusterings"]]) >=1)) {
      #to create this list of clusterings, use function "heatmap_artifacts_from_Seurat" below(output=markers, all_markers_lognorm, etc.)
      write_h5_clusterings(fn, artifact_name, layers[[artifact_name]][["clusterings"]])
    }
  }
}


#' Write one layer into an h5 file
#'
#' @description
#' Write an iscvam compliant layer into an h5 file.
#' @param fn This is the name of the h5 file.
#' @param artifact_name This is the name of the layer.
#' @param artifact This is a list of analytical entities for this dataset related to the specified layer.
#'     Required elements are:
#'         covs, a cell by covariate matrix of covariates.
#'             columns tsne_1, tsne_2, umap_1, umap_2 are required.
#'         discreteCovs, a vector of strings designating columns in covs that are discrete.
#'         continousCovs, a vector of strings designating columns in covs that are continuous.
#'
#' @examples
#'
#'\dontrun{
#'     write_h5_artifact("my.iscvam.h5", "myeloid", myeloid.artifact)
#' }
#' @export
write_h5_artifact <- function(fn, artifact_name, artifact) {
  # Check if file exists
  stopifnot(file.exists(fn))

  # Create a group in the HDF5 file
  rhdf5::h5createGroup(fn, paste0("artifacts/", artifact_name))

  # Convert factor covs to character
  is_factor <- sapply(artifact$covs, is.factor)
  artifact$covs[is_factor] <- lapply(artifact$covs[is_factor], as.character)

  # Write covs, discreteCovs, and continuousCovs to HDF5 file
  artifact_attributes <- c("covs", "discreteCovs", "continuousCovs")
  for (attribute in artifact_attributes) {
    rhdf5::h5write(artifact[[attribute]],
                   fn,
                   paste0("artifacts/", artifact_name, "/", attribute))
  }
}



#' Write a collection of clustering artifacts for a specific layer into an h5 file
#'
#' @description
#' Write a collection of iscvam compliant clustering artifacts for a specified layer into an h5 file.
#' @param fn This is the name of the h5 file.
#' @param layer This is the name of the layer.
#' @param clusterings This is a list of clustering artifacts related to the specified layer.
#'
#' @examples
#'\dontrun{
#'    write_h5_clusterings("my.iscvam.h5", "myeloid", myeloid.clusterings)
#'}
#' @export
write_h5_clusterings <- function(fn, layer, clusterings) {
  # Check if file exists
  stopifnot(file.exists(fn))

  # Create a group in the HDF5 file
  rhdf5::h5createGroup(fn, paste0("artifacts/", layer, '/clusterings'))

  # For each clustering scheme, save heatmap artifacts
  for (clustering in names(clusterings)) {
    path <- paste0("/artifacts/", layer, "/clusterings/", clustering)

    rhdf5::h5createGroup(fn, path)

    current_clustering <- clusterings[[clustering]]

    # Heatmap artifacts can be multimodal, in which case, one set of heatmap artifacts
    # may be stored for each modality.
    if ("markers" %in% names(current_clustering)) {
      write_h5_clustering(fn, path, current_clustering)
    } else {
      for (modality in names(current_clustering)) {
        rhdf5::h5createGroup(fn, paste0(path, "/", modality))
        write_h5_clustering(fn, paste0(path, "/", modality), current_clustering[[modality]])
      }
    }
  }
}


#' Write one set of clustering artifacts for a specific layer into an h5 file
#'
#' @description
#' Write a collection of iscvam compliant clustering artifacts for a specified layer into an h5 file.
#' @param fn This is the name of the h5 file.
#' @param path The location in the h5 file to store the heatmap artifacts related to a specific clustering scheme
#' @param clustering This is a list of heatmap artifacts related to the specified clustering scheme.
#'     Required elements:
#'         markers, a matrix of differential expression results
#'         all_markers_scaled, a matrix of per cluster average expression as z scores.
#'         all_markers_lognorm, a matrix of per cluster average expression as lognorm of the fraction.
#'         heatmap_markers_scaled, a matrix of per cluster average expression as z scores, restricted to
#'             markers to be used in the heatmap
#'         heatmap_markers_lognorm, a matrix of per cluster average expression as lognorm of the fraction,
#'             restricted to markers to be used in the heatmap
#'
#' @examples
#'\dontrun{
#'    write_h5_clustering("my.iscvam.h5", "/artifacts/myeloid/clusterings/louvain_0.6",
#'        myeloid.louvain_0.6.clustering)
#'}
#' @export
write_h5_clustering <- function(fn, path, clustering) {
  # Ensure the file exists
  stopifnot(file.exists(fn))

  # Create a group in the HDF5 file
  # rhdf5::h5createGroup(fn, path)

  # Define a vector of heatmap artifact names
  heatmap_artifacts <- c("markers","all_markers_scaled",
                         "all_markers_lognorm", "heatmap_markers_scaled",
                         "heatmap_markers_lognorm")

  # Write each heatmap artifact to the HDF5 file
  for (artifact in heatmap_artifacts) {
    artifact_path <- paste0(path, "/", artifact)
    rhdf5::h5write(clustering[[artifact]], fn, artifact_path)
  }
}



#' Write an extra underlying assay into an h5 file
#'
#' @description
#' Write raw count data of an assay into an h5 file. The saved count data will become queryable in ISCVAM.
#' @param fn This is the name of the h5 file.
#' @param counts This is the feature by cell sparse matrix of assay data to be saved.
#' @param assay Name of the assay.
#'
#' @examples
#'\dontrun{
#'    write_h5_extra_assay("my.iscvam.h5", seurat@assays$peaks@counts, "peaks")
#'}
#' @export

write_h5_extra_assay <- function (fn, counts, assay) {
  # Check if file exists
  stopifnot(file.exists(fn))

  # Create root path
  root <- paste0('assaysData/', assay)

  # Create HDF5 groups
  rhdf5::h5createGroup(fn, root)
  rhdf5::h5createGroup(fn, paste0(root, "/matrix"))

  # Write data into HDF5 file
  rhdf5::h5write(colnames(counts), fn, paste0(root, "/matrix/barcodes"))
  rhdf5::h5write(counts@i, fn, paste0(root, "/matrix/indices"))
  rhdf5::h5write(counts@p, fn, paste0(root, "/matrix/indptr"))
  rhdf5::h5write(counts@x, fn, paste0(root, "/matrix/data"))
  rhdf5::h5write(counts@Dim, fn, paste0(root, "/matrix/shape"))
  rhdf5::h5write(counts@Dimnames[[1]], fn, paste0(root, "/matrix/gene_names"))
}

#' Build layer artifacts from a seurat object
#'
#' @description
#' Construct layer artifacts from a seurat object.
#' @param seurat The seurat object.
#' @param qc_features Names of the qc features stored in the seurat object.
#' @param clustering_names Names of the clustering schemes for this layer.
#' @param umap Name of the umap reductions for this layer.
#' @param tsne Name of the tsne reductions for this layer.
#' @param extra_discrete_covs Additonal discrete covariate matrix. Should have the same number of rows as seurat.
#' @param extra_continuous_covs Additonal continuous covariate matrix. Should have the same number of rows as seurat.
#' @examples
#'\dontrun{
#'    layer_artifacts_from_seurat(
#'        seurat,
#'        qc_feature=
#'            c("nFeature_RNA", "nCount_RNA",
#'                "nFeature_ATAC","nCount_ATAC",
#'                "nFeature_peaks","nCount_peaks",
#'                "percent.mt","percent.rp",
#'                "nucleosome_signal","nucleosome_percentile",
#'                "TSS.enrichment","TSS.percentile" ),
#'        clustering_names=c("seurat_clusters_0.6"))
#'}
#' @export
layer_artifacts_from_seurat <- function(seurat, qc_features, umap = "umap", tsne = "tsne",
                                        clustering_names, extra_discrete_covs = extra_discrete_covs,
                                        extra_continuous_covs = NULL) {

  # Extract UMAP coordinates and set column names
  umap_cords <- seurat@reductions[[umap]]@cell.embeddings
  colnames(umap_cords) <- c("umap_1", "umap_2")

  # Extract t-SNE coordinates and set column names
  tsne_cords <- seurat@reductions[[tsne]]@cell.embeddings
  colnames(tsne_cords) <- c("tsne_1", "tsne_2")

  # Get quality control stats
  qc_stats <- seurat[[qc_features]]

  # Extract clusterings and convert to numeric
  clusterings <- sapply(clustering_names, function(cf) {
    (as.character(unlist(seurat[[cf]])))
    # as.numeric(as.character(unlist(seurat[[cf]])))
  })

  # Combine covariates including IDs, coordinates, and QC stats
  if(is.null(extra_continuous_covs)) {
  covs <- cbind(extra_discrete_covs, clusterings, 
                id = names(seurat$orig.ident), tsne_cords, umap_cords, qc_stats)
  } else{
    covs <- cbind(extra_discrete_covs, clusterings, extra_continuous_covs,
                  id = names(seurat$orig.ident), tsne_cords, umap_cords, qc_stats)
  }

  # Identify discrete and continuous covariates
  discrete_covs <- c(colnames(extra_discrete_covs), clustering_names)
  continuous_covs <- c(colnames(extra_continuous_covs), colnames(tsne_cords),
                       colnames(umap_cords), colnames(qc_stats))

  # Return list with covariates and their types
  return(list(covs = covs, discreteCovs = discrete_covs, continuousCovs = continuous_covs))
}


#' Build heatmap artifacts for a clustering scheme from a seurat object
#'
#' @description
#' Construct iscvam compliant heatmap artifacts from a seurat object.
#' @param seurat The seurat object.
#' @param clustering A vector of cluster assignments.
#' @param markers Differential markers for this clustering scheme.
#' @param assay Assay to use.
#' @param min_lfg minimum log fold change to filter markers for heatmap
#' @param max_p_adj maximum adjusted p value to filter markers for heatmap
#' @param limit maximum number of markers to use for each cluster for heatmap
#' @param mask_mt_rp whether to mask mitocondrial and ribosomal genes for heatmap
#' @importFrom rlang .data

#' @examples
#'\dontrun{
#'    heatmap_artifacts_from_seurat(
#'        seurat, curated_cell_types, markers_curate_cell_types, "RNA")
#'}
#' @export
heatmap_artifacts_from_seurat <- function(seurat, clustering, markers, assay, min_lfg=0.0,
                                          max_p_adj = 0.05, limit=10, mask_mt_rp=TRUE) {
  Seurat::DefaultAssay(seurat) <- assay
  if(is.null(seurat[[assay]]@data)) {
    seurat <- Seurat::NormalizeData(seurat)
  }
  if(is.null(seurat[[assay]]@scale.data)) {
    seurat <- Seurat::ScaleData(seurat)
  }
  Seurat::Idents(seurat) <- clustering

  cluster_averages <- Seurat::AverageExpression(seurat, return.seurat = TRUE)


  markers_use <- markers %>% dplyr::filter(
    .data$avg_log2FC > min_lfg &
      .data$p_val_adj < max_p_adj &
      !(mask_mt_rp & (grepl("^(RP[SL][0-9])|(Rp[sl][0-9])", .data$gene) |
                        grepl("^(MT-)|(mt-)", .data$gene)))
  ) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::top_n(limit, -.data$p_val_adj) %>%
    dplyr::top_n(limit, abs(.data$pct.1-.data$pct.2)) %>%
    dplyr::arrange(.data$cluster)

  markers_use <- unique(as.character(markers_use$gene))

  extract_data_frame <- function(data_source, features) {
    as.data.frame(
      cbind(
        feature.name = rownames(data_source),
       # as.matrix(cluster_averages[[cluster_averages@active.assay]]@data[features,])
        
        as.matrix(data_source[features,])
      )
    )
  }

  all_markers_lognorm <- extract_data_frame(cluster_averages[[cluster_averages@active.assay]]@data, rownames(cluster_averages[[cluster_averages@active.assay]]@data))
  all_markers_scaled <- extract_data_frame(cluster_averages[[cluster_averages@active.assay]]@scale.data, rownames(cluster_averages[[cluster_averages@active.assay]]@scale.data))
  heatmap_markers_lognorm <- extract_data_frame(cluster_averages[[cluster_averages@active.assay]]@data[markers_use,], markers_use)
  heatmap_markers_scaled <- extract_data_frame(cluster_averages[[cluster_averages@active.assay]]@scale.data[markers_use,], markers_use)

  list(
    markers = markers %>% dplyr::mutate(cluster = as.character(.data$cluster)),
    all_markers_lognorm = all_markers_lognorm,
    all_markers_scaled = all_markers_scaled,
    heatmap_markers_lognorm = heatmap_markers_lognorm,
    heatmap_markers_scaled = heatmap_markers_scaled
  )
}



