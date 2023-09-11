
###############
### this script is an example of how to write ISCVAM-h5 compliant
## inputs will be R objects, can be either analysis results from running multiome pipeline, or curation R objects
# output will be an ISCVAM-h5 compliant file
# you can deposit the output h5 on ISCVAM to visualize analysis results

####################################################
load("seurat.wnn.clusterings.GeneActivity.RData")
load("clus.avg.ALL.layers_h5.RData")
load("clus.avg_main.cat_h5.RData")
load("clus.avg_sub.cat_h5.RData")
load("clus.avg_default.ct_h5.RData")
load("clus.avg_wnn_0.6_h5.RData")

load("clus.avg_curated.ct_h5_update.RData")
load("covs.discrete_update.reported.anno.RData")


##############################
source("gene_annotation.R")
source("iscvamR_h5.R")


####
#create heatmap artifacts, preparing for h5 file
heatmap_artifacts = assemble_heatmap_artifacts(clus.avg.ALL.layers_h5,
                                               clus.avg_main.cat_h5 = clus.avg_main.cat_h5,
                                               clus.avg_sub.cat_h5,
                                               clus.avg_default.ct_h5,
                                               clus.avg_wnn_0.6_h5,
                                               clus.avg_curated.ct_h5_update)


#write h5 multiome 
write_mm_h5(seurat, covs.discrete, heatmap_artifacts, "int.validation.add.gene.anno_3.h5")


#############################
###check h5 result

#file name
fn <- "int.validation.add.gene.anno_3.h5"

library(rhdf5)

#read h5 file
h5 <- H5Fopen(fn) 

#check our h5 structure
h5ls(h5)

#example of how you can read a specific clustering layer
markers <- h5$'/artifacts/rna/clusterings/sub.categories/rna/markers'


#close our h5 file
h5closeAll()


