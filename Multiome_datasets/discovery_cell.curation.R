###
#intermediate steps for cell annotation, before writing to h5-ISCVAM compliant
#input will be results from analysis pipeline
#very specific just for this DISCOVERY dataset


load("seurat.wnn.clusterings.RData")
load("covs.RData")

default.ct <- seurat.wnn.clusterings$default.ct


load("covs.discrete.RData")
load("clus.avg.ALL.layers_h5.RData")

load("clus.avg_sub.categories_update.RData") #var name is "clus.avg_sub.categories_test"
load("cluster_averages_default.ct_h5.RData") #var name is 'cluster_averages_default_ct_h5'

load("clus.avg_wsnn_0.6_h5.RData") #var name is "clus.avg_wsnn_0.6"
load("curated.ct_wnn_for_h5.RData") #var name is 'curated.ct_wnn_test'




clus.avg_wnn_0.6_h5 <- clus.avg_wsnn_0.6_h5
clus.avg_wnn_0.6_h5$rna$markers$cluster <- as.character(clus.avg_wnn_0.6_h5$rna$markers$cluster)

clus.avg_wnn_0.6_h5$atac$markers$cluster <- as.character(clus.avg_wnn_0.6_h5$atac$markers$cluster)

save(clus.avg_wnn_0.6_h5, file = "clus.avg_wnn_0.6_h5.RData")


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


clus.avg_sub.cat_h5 <- markers_convert.clus.to.string(clus.avg_sub.categories_test)
clus.avg_default.ct_h5 <- markers_convert.clus.to.string(cluster_averages_default_ct_h5)
clus.avg_curated.ct_h5 <- markers_convert.clus.to.string(curated.ct_wnn_test)


save(clus.avg_sub.cat_h5, file = "clus.avg_sub.cat_h5.RData")
save(clus.avg_default.ct_h5, file = "clus.avg_default.ct_h5.RData")
save(clus.avg_curated.ct_h5, file = "clus.avg_curated.ct_h5.RData")


load("clus.avg_sub.cat_h5.RData")
load("clus.avg_default.ct_h5.RData")
load("clus.avg_wnn_0.6_h5.RData")
load("clus.avg_curated.ct_h5.RData")


#################
## changing curated


old.curated <- names(table(covs.discrete$curated.ct_wnn))
covs.discrete$curated.ct_wnn <- gsub("Cycling", "Prolif", covs.discrete$curated.ct_wnn)
covs.discrete$curated.ct_wnn <- gsub("FOXP3.*","FOXP3+Exh-8", covs.discrete$curated.ct_wnn)

#####
#for clus.avg_curated.ct_h5, change all colnames 


change.col.names <- function(df){
  old.col.names <- colnames(df)
  
  new.col.names <- gsub("Cycling", "Prolif", old.col.names)
  new.col.names <- gsub("FOXP3.*","FOXP3+Exh-8",new.col.names)
  colnames(df) <- new.col.names
  return(df)
}

change.curated_markers <- function(marker.df){
  
  marker.df$cluster <- gsub("Cycling", "Prolif", marker.df$cluster)
  marker.df$cluster <- gsub("FOXP3.*","FOXP3+Exh-8", marker.df$cluster)
  return(marker.df)
}

clus.avg_curated.ct_h5_update <- clus.avg_curated.ct_h5

for (modality in names(clus.avg_curated.ct_h5_update)){
  curr_modality <- clus.avg_curated.ct_h5_update[[modality]]
  for (df.name in names(curr_modality)){
    df <- curr_modality[[df.name]]
    clus.avg_curated.ct_h5_update[[modality]][[df.name]] <- change.col.names(df)
    print(colnames(clus.avg_curated.ct_h5_update[[modality]][[df.name]]))
    
    if (df.name=="markers"){
      clus.avg_curated.ct_h5_update[[modality]][["markers"]] <- change.curated_markers(df)
      print(table(clus.avg_curated.ct_h5_update[[modality]][["markers"]]$cluster))
    }
    
  }
}

save(clus.avg_curated.ct_h5_update, file = "clus.avg_curated.ct_h5_update.RData")
#############
#for NK in sub.cat


covs.discrete$sub.categories <- gsub("NK.*", "NK-like cells", covs.discrete$sub.categories)


df <- clus.avg_sub.cat_h5[['rna']][["all_markers_lognorm"]]

change.col.names_sub.cat <- function(df){
  new.col.names <- gsub("NK.*", "NK-like cells", colnames(df))
  colnames(df) <- new.col.names
  return(df)
}

change.clus.markers_sub.cat <- function(df){
  df$cluster <- gsub("NK.*", "NK-like cells", df$cluster)
  return(df)
}



clus.avg_sub.cat_h5_update <- clus.avg_sub.cat_h5

for (modality in names(clus.avg_sub.cat_h5_update)){
  curr_modality <- clus.avg_sub.cat_h5_update[[modality]]
  for (df.name in names(curr_modality)){
    df <- curr_modality[[df.name]]
    clus.avg_sub.cat_h5_update[[modality]][[df.name]] <- change.col.names_sub.cat(df)
    print(colnames(clus.avg_sub.cat_h5_update[[modality]][[df.name]]))
    
    if (df.name=="markers"){
      clus.avg_sub.cat_h5_update[[modality]][["markers"]] <- change.clus.markers_sub.cat(df)
      print(table(clus.avg_sub.cat_h5_update[[modality]][["markers"]]$cluster))
    }
    
  }
}

save(clus.avg_sub.cat_h5_update, file = "clus.avg_sub.cat_h5_update.RData")

save(covs.discrete, file = "covs.discrete_update.RData")
#############

load("covs.discrete_update.RData")

seurat.wnn.clusterings$curated_cell_types <- gsub("Cycling", "Prolif", seurat.wnn.clusterings$curated_cell_types)
seurat.wnn.clusterings$curated_cell_types <- gsub("FOXP3.*","FOXP3+Exh-8", seurat.wnn.clusterings$curated_cell_types)

colnames(seurat.wnn.clusterings@meta.data)[91] <- "reported.annotation"


###run gene annotation in seurat
source("/Volumes/Chen_Y_Ann/forThanh/ISCVAM_pipeline/multiome/gene_annotation.R")

seurat.wnn.clusterings<- calculate_gene_activity(seurat.wnn.clusterings)

save(seurat.wnn.clusterings, file = "seurat.wnn.clusterings_update.anno.RData")

#########
#add reported.annotation in covs.discrete

#load("covs.discrete_update.RData")
covs.discrete$reported.annotation <- seurat[["reported.annotation"]]$reported.annotation

save(covs.discrete, file = "covs.discrete_update.reported.anno.RData")

#########################################





