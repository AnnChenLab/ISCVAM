######
# this script is an example of how users can update their cell annotation, or adding meta data to covs file
# inputs will be 8 Seurat R objects from multiome analysis pipeline
# ouputs will be the updated of those R objects with cell annotation across different clusterings layers

#this script is specific for internal validation dataset

###########
library(Seurat)
library(Signac)

#load R objects
load("covs.RData")
load("seurat.wnn.clusterings.RData")

load("clus.avg.ALL.layers_h5.RData")
load("clus.avg_main.categories_h5.RData")
load("clus.avg_sub.categories_h5.RData")
load("clus.avg_default.ct_h5.RData")
load("clus.avg_wnn_0.6_h5.RData")
load("clus.avg_curated.ct_h5.RData")

#############
#create covs.discrete from pipeline
default.ct <- seurat.wnn.clusterings$default.ct

covs <- get.covs(seurat.wnn.clusterings, singler.preds, "Human", default.ct)
seurat <- seurat.wnn.clusterings


covs.discrete <- cbind(curated.ct_wnn= seurat.wnn.clusterings$curated.ct, 
                       main.categories = seurat.wnn.clusterings$main.categories, 
                       sub.categories = seurat.wnn.clusterings$sub.categories, 
                       covs)


########################################
#modify covs.discrete, adding patient info to covs.discrete

sample.type <- gsub(".*-", "", seurat.wnn.clusterings@meta.data$dataset)
samples <- seurat.wnn.clusterings@meta.data$dataset
patient <- gsub("-.*", "",  seurat.wnn.clusterings@meta.data$dataset)


covs.update <- cbind(samples = samples, 
                     sample.type = sample.type, 
                     patient =patient, 
                     covs.discrete)

colnames(covs.discrete)[2] <- "sample.type"

covs.discrete <- covs.update
save(covs.discrete, file = "covs.discrete.RData")

############################################
## this section is to change formating of our clustering layers to prepare for writing h5 file

####
#in 6 clustering layers, need to convert 'factor' to 'character' column
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


clus.avg_main.cat_h5 <- markers_convert.clus.to.string(clus.avg_main.cat_h5)
clus.avg_sub.cat_h5 <- markers_convert.clus.to.string(clus.avg_sub.cat_h5)
clus.avg_default.ct_h5 <- markers_convert.clus.to.string(clus.avg_default.ct_h5)
clus.avg_curated.ct_h5 <- markers_convert.clus.to.string(clus.avg_curated.ct_h5)


clus.avg_main.cat_h5[["atac"]][["markers"]][["cluster"]] <- as.character(clus.avg_main.cat_h5[["atac"]][["markers"]][["cluster"]])
clus.avg_main.cat_h5[["rna"]][["markers"]][["cluster"]] <- as.character(clus.avg_main.cat_h5[["rna"]][["markers"]][["cluster"]])

clus.avg_sub.cat_h5[["atac"]][["markers"]][["cluster"]] <- as.character(clus.avg_sub.cat_h5[["atac"]][["markers"]][["cluster"]])
clus.avg_sub.cat_h5[["rna"]][["markers"]][["cluster"]] <- as.character(clus.avg_sub.cat_h5[["rna"]][["markers"]][["cluster"]])


save(clus.avg_main.cat_h5, file = "clus.avg_main.cat_h5.RData")
save(clus.avg_sub.cat_h5, file = "clus.avg_sub.cat_h5.RData")
save(clus.avg_default.ct_h5, file = "clus.avg_default.ct_h5.RData")
save(clus.avg_curated.ct_h5, file = "clus.avg_curated.ct_h5.RData")


###############
load("seurat.wnn.clusterings.GeneActivity.RData")
load("clus.avg.ALL.layers_h5.RData")
load("clus.avg_main.cat_h5.RData")
load("clus.avg_sub.cat_h5.RData")
load("clus.avg_default.ct_h5.RData")
load("clus.avg_wnn_0.6_h5.RData")

load("clus.avg_curated.ct_h5_update.RData")
load("covs.discrete_update.anno.RData")


###############################
#add original 'reported.annotation' in the paper to our covs.discrete

covs.discrete$reported.annotation <- seurat$reported.annotation
save(covs.discrete, file = "covs.discrete_update.reported.anno.RData")

