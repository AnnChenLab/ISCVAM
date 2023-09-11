
########
#write ISCVAM-h5 compliant

load("seurat.wnn.clusterings.GeneActivity.RData") #this R object named "seurat" instead of "seurat.wnn.clusterings"
load("clus.avg_curated.ct_h5_update.RData")
load("clus.avg_sub.cat_h5_update.RData")
load("covs.discrete_update.reported.anno.RData")
load("clus.avg_default.ct_h5.RData")
load("clus.avg_wnn_0.6_h5.RData")
load("clus.avg.ALL.layers_h5.RData")

#using ISCVAM R package, adding gene annotation
#in pipeline on github, this "gene_annotation.R" script is called "updated_write_h5_scripts_using_iscvamR.R"
source("gene_annotation.R")

source("iscvamR_h5.R")

heatmap_artifacts = assemble_heatmap_artifacts(clus.avg.ALL.layers_h5,
                                               clus.avg_main.cat_h5=NULL,
                                               clus.avg_sub.cat_h5_update,
                                               clus.avg_default.ct_h5,
                                               clus.avg_wnn_0.6_h5,
                                               clus.avg_curated.ct_h5_update)


write_mm_h5(seurat, covs.discrete, heatmap_artifacts, "discovery.github.h5")

