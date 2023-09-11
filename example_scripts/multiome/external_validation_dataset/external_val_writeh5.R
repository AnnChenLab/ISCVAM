######
#Example script for PBMC dataset. Since we don't do cell annotation for this data, 
#inputs will be R objects resulting from multiome analysis pipeline
#final output will be a single ISCVAM-h5 compliant file 


#load R objects from analysis pipeline
library(Seurat)
library(Signac)

load("seurat.wnn.clusterings.RData")
load("covs.discrete.RData")
load("clus.avg.ALL.layers_h5.RData")
load("clus.avg_main.categories_h5.RData")
load("clus.avg_sub.categories_h5.RData")
load("clus.avg_default.ct_h5.RData")
load("clus.avg_wnn_0.6_h5.RData")
load("clus.avg_curated.ct_h5.RData")


#################
#modify R objects format to prepare for writing h5 file
markers_convert.clus.to.string <- function(clus.avg_list){
  ### in marker df, convert cluster column from factor to character
  #so when writing h5 file, these cluster names will be conserved (instead of being converted to integer by h5)
  
  clus.avg_update <- lapply(names(clus.avg_list), function(modality){
    marker.df <- clus.avg_list[[modality]][["markers"]]
    marker.df$cluster <- as.character(marker.df$cluster)
    clus.avg_list[[modality]][["markers"]] <- marker.df
    return(clus.avg_list[[modality]])
  })
  names(clus.avg_update) <- names(clus.avg_list)
  return(clus.avg_update)
}


clus.avg_main.categories_h5 <- markers_convert.clus.to.string(clus.avg_main.categories_h5)
clus.avg_sub.categories_h5 <- markers_convert.clus.to.string(clus.avg_sub.categories_h5)
clus.avg_default.ct_h5 <- markers_convert.clus.to.string(clus.avg_default.ct_h5)
clus.avg_curated.ct_h5 <- markers_convert.clus.to.string(clus.avg_curated.ct_h5)

save(clus.avg_main.categories_h5, file = "clus.avg_main.categories_h5.RData")
save(clus.avg_sub.categories_h5, file = "clus.avg_sub.categories_h5.RData")
save(clus.avg_default.ct_h5, file = "clus.avg_default.ct_h5.RData")
save(clus.avg_curated.ct_h5, file = "clus.avg_curated.ct_h5.RData")

##########

load("seurat.wnn.clusterings_gene.annotation.RData")
load("clus.avg.ALL.layers_h5.RData")
load("clus.avg_main.categories_h5.RData")
load("clus.avg_sub.categories_h5.RData")
load("clus.avg_default.ct_h5.RData")
load("clus.avg_curated.ct_h5.RData")
load("clus.avg_wnn_0.6_h5.RData")
load("covs.discrete.RData")

#########################
#load R scripts to write h5 file
source("gene_annotation.R")
source("iscvamR_h5.R")



heatmap_artifacts = assemble_heatmap_artifacts(clus.avg.ALL.layers_h5,
                                               clus.avg_main.cat_h5 = clus.avg_main.categories_h5,
                                               clus.avg_sub.categories_h5,
                                               clus.avg_default.ct_h5,
                                               clus.avg_wnn_0.6_h5,
                                               clus.avg_curated.ct_h5)

##########################
#calculate gene activity (this step was not included in the analysis pipeline)
seurat.wnn.clusterings<- calculate_gene_activity(seurat.wnn.clusterings)

#########################
#write ISCVAM-h5 file 

write_mm_h5(seurat.wnn.clusterings, covs.discrete, heatmap_artifacts, "PBMC_multiome_github.h5")



#########################
#check the h5 result
library(rhdf5)


#file name
fn <- "PBMC_multiome_github.h5"

#read h5 structure
h5 <- H5Fopen(fn)
h5ls(h5)

#######################
#examples of checking h5 components

#checking covs variables
covs_h5 <- h5$'/artifacts/atac/covs'
colnames(covs_h5)

#check marker df
markers <- h5$'/artifacts/wnn/clusterings/curated.ct_wnn/rna/markers'
head(markers)


#close h5 file
h5closeAll()