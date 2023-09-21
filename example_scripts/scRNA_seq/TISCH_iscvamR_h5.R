############
#this scritp is to convert a Tisch dataset into ISCVAM-h5 compliant
#input will be Tisch's h5 file and "CellMetainfo_table"
#output will be ISCVAM-h5 file 

library(Seurat)
source("iscvamR_h5.R")


project <- "SKCM_GSE115978_aPD1"


###################
#create Seurat object and covs object, combining them into a "seurat.list" object
create.seurat.list <- function(data.folder, project){
  
  data.path <- paste0(data.folder, '/', project, "_Expression")
  
  #check files in folder
  #list.files(data.path)
  
  
  #read raw 10x file and create Seurat object
  data <- Read10X_h5(paste0(data.path, '/', project, "_expression.h5"))
  seurat.raw <- CreateSeuratObject(counts=data,project=project)
  
  #create meta.data file
  meta.file <- list.files(data.path, pattern =".*_CellMetainfo_table.tsv")
  meta.data <- read.table(file = paste0(data.path, '/', meta.file), 
                          sep = '\t', header = TRUE)
  
  #create Seurat list
  seurat.list <- list(all=list(seurat = seurat.raw, 
                               covs = meta.data))
  return(seurat.list)
}

seurat.list <- create.seurat.list(data.folder, project)


######################
#add cell annotation to Seurat object
seurat <- seurat.list$all$seurat
covs <- seurat.list$all$covs


add.cell.anno.to.seurat.Tisch <- function(seurat, covs){
  ## add cell.annotation (Celltype) and Cluster numbers (Cluster) to seurat from covs
  clusterings <- grep("(Celltype.*)|(Cluster)", colnames(covs), value=TRUE)
  clusterings.df <- covs[clusterings]
  rownames(clusterings.df) <- covs$Cell
  
  #add clusterings to seurat
  seurat <- AddMetaData(seurat, clusterings.df)
  return(seurat)
}

seurat <- add.cell.anno.to.seurat.Tisch(seurat, covs)

###########################
##find markers for cell annotation in Tisch 
find.markers.by.all.clusterings_Tisch <- function(seurat, clusterings) {
  
  ret <- lapply(clusterings, function(clustering) {
    find.markers.by.clustering(seurat,  seurat[[clustering]])
  })
  names(ret) <- clusterings
  ret
}

all.markers <- find.markers.by.all.clusterings_Tisch(seurat, clusterings) 

######################
#combine all clusterings and their markers into "artifacts" for ISCVAM

create.artifacts.ALL.clusterings <- function(seurat, all.markers){
  ##create artifacts for ALL clusterings
  clusterings_artifacts <- lapply(names(all.markers), function(clustering){
    print(clustering)
    artifacts <- heatmap_artifacts_from_seurat(seurat=seurat, 
                                               clustering = clustering, 
                                               markers = all.markers[[clustering]], 
                                               assay="RNA")
  })
  names(clusterings_artifacts) <- names(all.markers)
  return(clusterings_artifacts)
}

clusterings_artifacts <- create.artifacts.ALL.clusterings(seurat, all.markers)


#########################
#customized Tisch covs for ISCVAM

prep_covs_Tisch <- function(covs){
  #this function is to (1) strip all non-clinical covs and only keep meta.data from cohort
  #(2) prep for input of artifacts in covs to write h5
  
  #remove cell barcodes in first column
 covs <- covs[,-1]

  qc_features = c("nFeature_RNA", "nCount_RNA")
  clustering_names = grep("(Cluster)|(Celltype.*)",colnames(covs), value=T)
  
  
  covs.clinical <- covs[, !colnames(covs) %in% c(qc_features, clustering_names)]
  
  return(list(qc_features = qc_features, 
              clustering_names = clustering_names,
              extra_discrete_covs = covs.clinical))
}

covs.Tisch <- prep_covs_Tisch(covs)

#############################
#create RNA layer artifacts
layer_artifacts_from_seurat_TISCH <- function(seurat, qc_features,
                                              clustering_names, extra_discrete_covs = NULL,
                                              extra_continuous_covs = NULL) {
  #TISCH doesnt have TSNE projection
  
  # Extract UMAP coordinates and set column names
  umap_cols = grep("UMAP_.*",colnames(extra_discrete_covs), value=T)
  umap_cords = extra_discrete_covs[umap_cols]
  colnames(umap_cords) <- c("umap_1", "umap_2")
  
  #then remove umap_cols from extra_discrete_covs
  extra_discrete_covs <- extra_discrete_covs[, !colnames(extra_discrete_covs) %in% umap_cols]
  
  # Get quality control stats
  qc_stats <- seurat[[qc_features]]
  
  # Extract clusterings and convert to numeric
  clusterings <- sapply(clustering_names, function(cf) {
    as.character(unlist(seurat[[cf]]))
    # as.numeric(as.character(unlist(seurat[[cf]])))
  })
  
  clusterings <- data.frame(clusterings)
  
  # Combine covariates including IDs, coordinates, and QC stats
  
  if(is.null(extra_continuous_covs)) {
    covs.h5 <- cbind(extra_discrete_covs, clusterings, 
                     id = colnames(seurat), umap_cords, qc_stats)
  }else{
    covs.h5 <- cbind(extra_discrete_covs, clusterings, extra_continuous_covs,
                     id = colnames(seurat), umap_cords, qc_stats )
  }
  
  # Identify discrete and continuous covariates
  discrete_covs <- c(colnames(extra_discrete_covs), clustering_names)
  continuous_covs <- c(colnames(extra_continuous_covs),
                       colnames(umap_cords), colnames(qc_stats))
  
  # Return list with covariates and their types
  list(covs = covs.h5, discreteCovs = discrete_covs, continuousCovs = continuous_covs)
}

covs_artifacts <- layer_artifacts_from_seurat_TISCH(seurat, 
                                                    qc_features = covs.Tisch$qc_features, 
                                                    clustering_names = covs.Tisch$clustering_names, 
                                                    extra_discrete_covs = covs.Tisch$extra_discrete_covs)

#adding clustering artifacts into covs
covs_artifacts[["clusterings"]] <- clusterings_artifacts


layers <- list(all = covs_artifacts)

################################
###writing ISCVAM-h5 file

fn <- paste0(project, "_test_6.h5")

write_h5(fn, seurat, layers, assays = c("RNA"))

############################
#check ISCVAM-h5 result
library(rhdf5)

h5 <- H5Fopen(fn)

#check h5 structure 
h5ls(h5)

#check variable names in covs
colnames(h5$'/artifacts/all/covs')

#check discrete covs 
h5$'/artifacts/all/discreteCovs'

