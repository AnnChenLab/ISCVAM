library(rhdf5)
library(Seurat)
library(dplyr)
library(SingleR)
library(future)


##################################


get.sample.table <- function(config) {
  # to be edited for sample information
  rna_data <- list.files(config$data, pattern = ".*RNA-",recursive = F)
  
  samples <- sapply(strsplit(rna_data,"-"), function(s){paste(head(s,-1), collapse = "-")})
 
  count.dirs.keep <- paste0(rna_data, "/outs/filtered_feature_bc_matrix.h5")
  sample.table <- data.frame(sample=samples, path= paste0(config$data, '/', count.dirs.keep) , stringsAsFactors = F)
  write.table(sample.table, paste0(config$stem,".sample.table"))
  #remove 1st line, 10-E-RNA since data not avai
  sample.table <- sample.table[-1,]
  sample.table
}


#sample.table <- get.sample.table(config)



################################
load.reads <- function(path.10x) {
  require(Seurat)
  # working only with 10x
  if(endsWith(path.10x, ".h5")) {
    scData = Read10X_h5(path.10x)
    if('Gene Expression' %in% names(scData)) {
      # seurat version changes
      scData <- scData[['Gene Expression']]
    }
  } else {
    scData = Read10X(path.10x)
  }
  scData
}

#############################################
#read h5 from batch ######################

load.reads.batch <- function(sample.table) {
  print(sample.table$path)
  ret <- lapply(sample.table$path, load.reads)
  #names(ret) <- sample.table$samples
  ret <- setNames(ret, sample.table$sample)
  ret
}


#raw.reads.batch <- load.reads.batch(sample.table)


######################################
##### create Seurat object #######

get.seurat.from.batch <- function(reads.batch, stem) {
  require(Seurat)
  s <- lapply(names(reads.batch), function(rb) {
    ret <- CreateSeuratObject(counts=reads.batch[[rb]],project=stem)
    ret$sample=rb
    ret
  })
  s.merge <- merge(s[[1]], y=tail(s,-1), add.cell.ids = 1:length(s), project=stem)
}


#seurat.raw <- get.seurat.from.batch(raw.reads.batch, config$stem)
#############################
###### create seurat.raw from TISCH 

create.seurat.list.from.TISCH <- function(data.path, dataset){
  
  #read raw 10x file and create Seurat object
  data <- Read10X_h5(paste0(data.path, '/', dataset, "_expression.h5"))
  seurat.raw <- CreateSeuratObject(counts=data,dataset=dataset)
  
  #create meta.data file
  meta.file <- list.files(data.path, pattern =".*_CellMetainfo_table.tsv")
  meta.data <- read.table(file = paste0(data.path, '/', meta.file), 
                          sep = '\t', header = TRUE)
  
  #create Seurat list
  seurat.list <- list(all=list(seurat = seurat.raw, 
                               covs = meta.data))
  return(seurat.list)
}

get.seurat.raw <- function(seurat.raw.list.from.TISCH){
  seurat.raw <- seurat.raw.list.from.TISCH$all$seurat
  return(seurat.raw)
}
###################

do.singler <- function(raw.reads, single.ref.rds) {
  require(SingleR)
  require(future)
  load(single.ref.rds)
  options(future.globals.maxSize = 30000 * 1024^2)
  preds <- sapply(names(refs), function(ref) {
    sapply(c("main","fine"), function(r) {
      SingleR(raw.reads, refs[[ref]], refs[[ref]][[paste0("label.",r)]] ) }
    )
  })
}

#singler.preds <- do.singler(seurat.raw@assays$RNA@counts, config$single.ref.rds)


#############################

qc.seurat <- function(seurat, species, nFeature=500) {
  mt.pattern <- case_when(
    species =="Human" ~ "^MT-",
    species =="Mouse" ~ "^mt-",
    TRUE ~ "^MT-"
  )
  ribo.pattern <- case_when(
    species == "Human" ~ "^RP[LS]",
    species == "Mouse" ~ "^Rp[ls]",
    TRUE ~ "^RP[LS]"
  )
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = mt.pattern)
  seurat[["percent.rp"]] <- PercentageFeatureSet(seurat, pattern = ribo.pattern)
  seurat[,seurat[["percent.mt"]] <= 20 & seurat[["nFeature_RNA"]] >= nFeature & seurat[['nCount_RNA']]>=1000]
  
}


  

###############################
#debatch
integrate.seurat <- function(seurat.qced, nGenes, nCCA){ 
  # nCCA=30
  # nGenes=6000

  read_depth <- ifelse(seurat.qced$sample %in% c( "4-E-RNA","7-E-RNA", "10-E-RNA", "14-E-RNA", "14-1-RNA"), "high_read_depth", "low_read_depth")
  seurat.qced$read_depth <- read_depth

  samples <- SplitObject(seurat.qced, split.by = "read_depth") %>% lapply(SCTransform)
  features <- SelectIntegrationFeatures(object.list = samples, 
                                        anchorfvf.nfeatures = 4000, # nfeatures for FindVariableFeatures. Used if VariableFeatures have not been set for any object in object.list.
                                        nfeatures = nGenes)

  samples <- PrepSCTIntegration(object.list = samples, anchor.features = features)

  require(future)
  options(future.globals.maxSize = 65000 * 1024^2)
  plan("multiprocess", workers = 6)

  anchors <- FindIntegrationAnchors(object.list = samples, anchor.features = features,
                                           normalization.method = "SCT", dims = 1:nCCA)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  print('done with integrate data')
  return(integrated)
}



analyze.seurat <- function(seurat, pcs.to.compute=200, pcs.to.analyze=40, debatch=FALSE) {
  require(Seurat)
  require(dplyr)
 
  #skip this step if already DEBATCH, since SCT already run in debatching steps
  if(debatch==FALSE){
    DefaultAssay(seurat) <- "RNA"
    seurat <- seurat %>% SCTransform(vars.to.regress = c("percent.mt", "percent.rp"),conserve.memory = T)
  }
 
  seurat <- seurat %>% RunPCA(npcs=pcs.to.compute)
  seurat <- seurat %>% FindNeighbors(dims=1:pcs.to.analyze)
  seurat <- seurat %>% FindClusters()
  seurat <- seurat %>% RunUMAP(dims = 1:pcs.to.analyze)
  seurat <- seurat %>% RunTSNE(dims = 1:pcs.to.analyze)
}



do.clusterings <- function(seurat, pcs.to.analyze) {
  for(resolution in c(0.6,0.8, 1, 1.2, 2, 4)) {
    seurat <- FindClusters(seurat, resolution = resolution)
    seurat[[paste0("seurat_clusters_",resolution)]] = seurat$seurat_clusters
  }
  return(seurat)
}

#####################################
#########read results of singleR #####

subset.preds <- function (preds, cells)  {
  #this function was used to create an EMPTY df to store singleR result'''
  levels <- rownames(preds)
  names <- colnames(preds)
  ret <- array(list(), c(length(levels), length(names)))
  rownames(ret) = levels
  colnames(ret) = names
  for(refn in names) {
    for(refl in levels) {
      ret[[refl, refn]] = preds[[refl,refn]][cells,]
    }
  }
  ret
}

preds.to.labels <- function(preds) {
  
  #this function is to pass the result of singleR to empty df
  
  ref.levels <- c("main", "fine")
  ref.names <- colnames(preds)
  ret <- list()
  for(refn in ref.names) {
    for(refl in ref.levels) {
      ret[[paste0(refn,".", refl)]] = preds[[refl,refn]]$pruned.labels
    }
  }
  ret
}


########################
#create covs

# get.covs.clinical<- function(seurat.analyzed){
#   #must change according to projects
#      
#     sample <- seurat.analyzed@meta.data$sample
#     covs <- as.data.frame(sample)
#     
#     #covs$sample <- gsub(-1-RNA, "-I-RNA", covs$sample)
# 
#     covs$patient <- gsub("-[E|I]-RNA", "", covs$sample)
# 
# 
#     Date_of_collection <- case_when(
#       grepl("4-I-RNA", covs$sample) ~ "9/27/19", 
#       grepl("4-E-RNA", covs$sample) ~ "2/13/20", 
#       grepl("5-I-RNA", covs$sample) ~ "2/19/20", 
#       grepl("5-E-RNA", covs$sample) ~ "10/8/20",
#       grepl("7-I-RNA", covs$sample) ~ "4/29/20", 
#       grepl("7-E-RNA", covs$sample) ~ "5/21/20", 
#       grepl("10-E-RNA", covs$sample) ~ "11/20/20", 
#       grepl("14-I-RNA", covs$sample) ~ "4/22/21", 
#       grepl("14-E-RNA", covs$sample) ~ "5/18/21"
#     )
# 
#     response <- case_when(
#       grepl("^4$", covs$patient) ~ "Responder", 
#       grepl("5", covs$patient) ~ "Responder", 
#       grepl("7", covs$patient) ~ "Non-responder", 
#       grepl("10", covs$patient) ~ "Non-responder",
#       grepl("^14$", covs$patient) ~ "Non-responder"
#     )
# 
#     disease <- case_when(
#       grepl("^14$", covs$patient) ~ "Pancreatic Cancer w/ LMD", 
#       TRUE ~ "Breast Cancer w/ LMD"
#     )
# 
#     mutation <- case_when(
#       grepl("^4$|7|^14$", covs$patient) ~ "ER/PR+, Her2-",
#       grepl("^5$|^10$", covs$patient) ~ "Triple Negative"
#     )
# 
#     pre_post <- case_when(
#       grepl(".*-I-.*", covs$sample) ~ "pre-treatment", 
#       grepl(".*-E-.*", covs$sample) ~ "post-treatment"
#     )
# 
# 
#     covs <- cbind(covs, Date_of_collection, response, disease, mutation, pre_post)
# 
# }

get.covs <- function(seurat, preds, config, covs.clinical) {
  ref.by.species=list(Mouse=c("immgen", "mousernaseq"), Human=c("blueprint","hpca", "dice", "monaco","novershtern"))
  #this function is to get 10 cols of cell type prediction from singleR
  preds.labels <- preds.to.labels(
                          subset.preds(preds,colnames(seurat))[,ref.by.species[[config$species]]])

  covs <- cbind(covs.clinical, preds.labels)
  return(covs)
}

get.covs_debatch <- function(seurat, covs,  annotation, layer) {
  #here dont use preds.labels, because this is DEBATCH, colnames(seurat) of debatch already merge and mess up, 
  #thus when use colnames(seurat) to retrieve singleR, will not get the right values
  
  # ref.by.species=list(Mouse=c("immgen", "mousernaseq"), Human=c("blueprint","hpca", "dice", "monaco","novershtern"))
  # preds.labels <- preds.to.labels(
  #   subset.preds(preds,colnames(seurat))[,ref.by.species[[config$species]]])
  
  covs_layer <- covs[covs$main.categories %in% layer,]
  ref_panel <- covs_layer %>% select(contains(c("blueprint","hpca", "dice", "monaco","novershtern")))
  covs = data.frame(
    sample=seurat$sample,
    curated.cell.types=annotation$curated.cell.types,
    ref_panel
  )
  covs
}

get.covs_2ndstage<- function(covs, layer, annotation){
  
  covs_layer <- covs[covs$main.categories %in% layer,]
  covs_layer <- subset(covs_layer, select = -c(main.categories,sub.categories,curated.cell.types))
  require("tibble")
  covs_layer <- add_column(covs_layer, curated.cell.types=annotation$curated.cell.types, .after = "read_depth")
  
}
###############

find.markers.by.clustering <- function(seurat, clusters) {
  DefaultAssay(seurat) <- "RNA"
  Idents(seurat) <- clusters
  
  #enable parallelization
  require(future)
  options(future.globals.maxSize = 64000 * 1024^2)
  plan("multisession", workers = 4)
  
  markers <- FindAllMarkers(seurat)
  markers$diff <- markers$pct.1 - markers$pct.2
  return(markers)
}


find.markers.by.all.clusterings <- function(seurat) {
  clusterings <- c("seurat_clusters_0.6",
                  "seurat_clusters_0.8", 
                   "seurat_clusters_1", 
                   "seurat_clusters_1.2",
                   "seurat_clusters_2",
                   "seurat_clusters_4")
  ret <- lapply(clusterings, function(clustering) {
    find.markers.by.clustering(seurat,  seurat[[clustering]])
  })
  names(ret) <- clusterings
  ret
}

#####################################
###write all marker to xlsx 

write_all.markers <- function(fileName, all.markers){
  library(openxlsx)
  
  wb <- createWorkbook(fileName)
  
  for (cluster in names(all.markers)){
    sheetName = cluster
    markers = all.markers[[cluster]]
    # Add worksheets to workbook
    addWorksheet(wb, sheetName = sheetName)
    writeData(wb, sheet = sheetName, markers)
  }
  
  saveWorkbook(wb, file = fileName, overwrite = TRUE)
}
####################


get.cell.and.cluster.stats.for.clustering <- function(seurat, preds, clustering="seurat_clusters") {
  #this function only use blueprint.main col 
  ct.main <- subset.preds(preds,colnames(seurat))[["main", "blueprint"]]$pruned.labels
  
  ncell.by.ct.main <- table(ct.main)
  ncell.by.ct.main <- sort(ncell.by.ct.main, decreasing = T)
  top.cts <- names(ncell.by.ct.main[ncell.by.ct.main>0])
  cell.stats <- cbind(seurat@meta.data, ct.main)
  
  cell.stats$clusterN = unlist(seurat[[clustering]])
  cell.stats$ct.main = factor(ct.main, levels=top.cts)
  
  cluster.stats <- cell.stats %>% 
    group_by(clusterN, ct.main) %>%
    summarize(n = n(), nFeature_RNA=median(nFeature_RNA), nCount_RNA=median(nCount_RNA),
              percent.mt = median(percent.mt)) %>%
    # count(clusterN, Blueprint ) %>%
    group_by(clusterN) %>%
    mutate(freq=round(n/sum(n), digits = 2))
  
  majority.stats <- cluster.stats %>% group_by(clusterN) %>% filter(freq==max(freq))
  list(cell.stats=cell.stats, cluster.stats=majority.stats)
}

get.cell.and.cluster.stats.for.all.clusterings <- function(seurat, preds) {
  clusterings <- c("seurat_clusters_0.6",
                    "seurat_clusters_0.8", 
                   "seurat_clusters_1", 
                   "seurat_clusters_1.2",
                   "seurat_clusters_2",
                   "seurat_clusters_4")
                  # "clusters_infomap")
  sapply(clusterings, function(clustering) {
    get.cell.and.cluster.stats.for.clustering(seurat, preds,clustering)
  })
}


write_cells.clusters.stats <- function(fileName, cell.and.cluster.stats){
  require(openxlsx)
  excel <- createWorkbook(fileName)
  
  #write multiple tabs/worksheets
  clusterings <- c("seurat_clusters_0.6",
                   "seurat_clusters_0.8", 
                   "seurat_clusters_1", 
                   "seurat_clusters_1.2",
                   "seurat_clusters_2",
                   "seurat_clusters_4")
  
  sapply(clusterings, function(clustering) {
    cell.stats <- cell.and.cluster.stats[["cell.stats", clustering]]
    cluster.stats <- cell.and.cluster.stats[["cluster.stats",clustering]]
    cluster_name <- gsub('seurat_', '', clustering)
    
    # Add worksheets to workbook
    addWorksheet(excel, sheetName = paste0(cluster_name, '_cell.stats'))
    writeData(excel, sheet = paste0(cluster_name, '_cell.stats'), cell.stats)
    
    addWorksheet(excel, sheetName = paste0(cluster_name, '_cluster.stats'))
    writeData(excel, sheet = paste0(cluster_name, '_cluster.stats'), cluster.stats)
  })
  
  saveWorkbook(excel, file = fileName, overwrite = TRUE)
}

#################################################
do.curated.cell.types <- function(cell.and.cluster.stats) {
  clustering = "seurat_clusters_0.6"
  cell.stats <- cell.and.cluster.stats[["cell.stats", clustering]]
  cluster.stats <- cell.and.cluster.stats[["cluster.stats", clustering]]
  
  major.cluster.cell.types <- cluster.stats$ct.main
  names(major.cluster.cell.types) <- cluster.stats$clusterN
  majority.exception.clusters = c()
  cluster <- cell.stats[[clustering]]
  curated.cell.types <- case_when(
    (! (cluster %in% majority.exception.clusters)) ~
      paste0(as.character(major.cluster.cell.types[as.character(cluster)]), "-", cluster),
    TRUE ~ "UNKNOWN"
  )
  curated.cell.types
}
####################################

annotate.curated.cell.types <- function(curated.cell.types, markers) {
  cts = unique(curated.cell.types)
  #extract last hyphen and cluster numbers of cts
  library(dplyr)
  library(stringr)
  clus_numb <- str_extract(cts, "-[0-9]+$")
  clus_numb <- gsub("-", "", clus_numb)
  
  ct <- as.character(strsplit(cts, "-[0-9]+$"))
  
  df <- data.frame(cts, clus_numb, ct) 
  df$clus_numb <- as.numeric(df$clus_numb)
  #reorder cts
  df <- df %>% arrange(clus_numb)
  cts.order <- df$cts
  
  markers$cluster <- factor(markers$cluster, levels = cts.order)
  markers <- markers %>% arrange(cluster)
  
  return(list(curated.cell.types=curated.cell.types, markers=markers))
}



###########################
plot.heatmap.seurat <- function(seurat, markers, fileName) {
  library(Seurat)
  pdf(file= fileName, height=28,width=20)
  
  DefaultAssay(seurat) <- "RNA"
  
  seurat <- NormalizeData(seurat, assay = "RNA")
  
  seurat <- ScaleData(seurat)
  
  cluster.averages <- AverageExpression(seurat, return.seurat = TRUE)
  
  require(ggplot2)
  require(dplyr)
  min.lfg <- 0.0
  max.p.adj <- 0.05
  
  n.clusters <- length(unique(markers$cluster))
  limit <- 140/n.clusters
  
  for(mask.mt.rp in c(T,F)) {
    markers.use=subset(markers,avg_log2FC > min.lfg & 
                         p_val_adj < max.p.adj & 
                         !(mask.mt.rp & 
                          (grepl("^RP[SL]", gene) | grepl("^MT", gene)))) %>%
      group_by(cluster)%>% 
      top_n(limit, -p_val_adj) %>% 
      top_n(limit, abs(pct.1-pct.2)) %>%
      
      
      #if clusters is in FACTOR, then "arrange(as.factor(cluster))"
      #because as.character can't sort alpha_numeric (cant sort 10 instead of 1)
      arrange(as.factor(cluster))
     # arrange(as.character(cluster))
    markers.use <- as.character(markers.use$gene)
    
    print(DoHeatmap(cluster.averages, features = markers.use,
                    draw.lines = FALSE) + 
            ggtitle(paste0("Average expression; mask.mt.rp=", mask.mt.rp))) 
         
    
    print(DoHeatmap(seurat, features = markers.use)+
            ggtitle(paste0("test=wilcox; zscore; mask.mt.rp=", mask.mt.rp)))
          
    
    print(DoHeatmap(seurat, features = markers.use, slot="data" )+
            ggtitle(paste0("test=wilcox; log normalized; mask.mt.rp=", mask.mt.rp)))
  
    
    if(nrow(seurat@assays[[seurat@active.assay]]@counts) > 0) {
      print(DoHeatmap(seurat, features = markers.use, slot="counts" )+
              ggtitle(paste0("test=wilcox; counts; mask.mt.rp=", mask.mt.rp)))
    }
  }
  dev.off()
}

#######################
create.report <- function(seurat, preds, config, version='original', marker.sets) {
  pdf(file=paste0(version, ".report.pdf"), height=30,width=24)
  print(seurat %>% ElbowPlot(ndims=config$pcs.to.compute))
 
  curated.cell.types = marker.sets$curated.cell.types
  Idents(seurat) = curated.cell.types
  levels(seurat) <- names(table(marker.sets$markers$cluster))
  
  print(seurat %>% VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 2))
  plot1 <- seurat %>% FeatureScatter(feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- seurat %>% FeatureScatter(feature1 = "nCount_RNA", feature2 = "percent.rp")
  plot3 <- seurat %>% FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(CombinePlots(plots = list(plot1, plot2, plot3)))
  print(DimPlot(seurat, reduction = "umap"))
  print(DimPlot(seurat, reduction = "tsne"))
  print(PCAPlot( seurat))
   
  dev.off()

  fileName <- paste0(version, ".report.heatmap.","seurat_0.6", ".pdf")
  plot.heatmap.seurat(seurat, marker.sets$markers, fileName)
    
}

##########################

create.seurat.list.for.ISCVAM <- function(seurat.analyzed, covs){
          seurat.list <- list(all=list(seurat = seurat.analyzed, 
                                       covs = covs))
          return(seurat.list)
          }

#########################
do.main.categories.lymphoids <- function(cts) {
  case_when(
   
    grepl("^NK cells", cts) ~ "NK cells",
    grepl("^B-cells", cts) ~ "B-cells",
    grepl("^T-cells", cts) ~ "T-cells",
    grepl("^CD4", cts) ~ "CD4",
    grepl("^CD8", cts) ~ "CD8",
    grepl("^pDC", cts) ~ "pDC", 

    TRUE ~ cts
  )
}

do.main.categories.tumor <- function(cts){
  case_when(
    grepl("^Epithelial cells", cts) ~ "Epithelial cells",
    grepl("^Neurons", cts) ~ "Neurons"
  )
}




#########################
do.subset.analyses <- function(seurat,  cells,  pcs.to.compute=200, pcs.to.analyze=40) {
  analyze.seurat(seurat[,cells],  pcs.to.compute=200, pcs.to.analyze=40)
}


