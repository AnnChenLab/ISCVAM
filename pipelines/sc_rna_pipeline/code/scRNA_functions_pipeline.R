library(rhdf5)
library(Seurat)
library(dplyr)
library(SingleR)
# library(future)
library(janitor)


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

###########
#create seurat object from Tisch 10x h5 
create.seurat.from.Tisch.h5 <- function(data.path, dataset.1){
  print(paste0("Creating Seurat object from dataset: ", dataset.1))
  data.folder <- paste0(data.path, '/', dataset.1)
  list.files(data.folder)
  h5.name <- list.files(data.folder, pattern = '.h5')
  counts.matrix <- Read10X_h5(paste0(data.folder, '/', h5.name), use.names = TRUE, unique.features = TRUE)
  
  seurat <- CreateSeuratObject(counts.matrix, project=dataset.1)
  seurat@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]] <- counts.matrix@Dimnames[[1]]
  seurat@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]] <- counts.matrix@Dimnames[[2]]
  
  meta.fn <- list.files(data.folder, pattern='_CellMetainfo_table.tsv')
  meta.data <- read.delim(paste0(data.folder, '/', meta.fn))
  
  #add meta data to seurat
  seurat@meta.data <- cbind(seurat@meta.data, meta.data)
  
  print(paste0("dataset dimension before qc, number of genes: ", dim(seurat)[1]))
  print(paste0("dataset dimension before qc, number of cells: ", dim(seurat)[2]))
  return(seurat)
}

#########
remove_umap_from_tisch <- function(seurat.raw){
  #remove pre-computed UMAP coord in tisch
  seurat.raw@meta.data <- seurat.raw@meta.data[,!colnames(seurat.raw@meta.data) %in% c("UMAP_1","UMAP_2")]
  return(seurat.raw)
}

###################

do.singler <- function(raw.reads, single.ref.rds) {
  require(SingleR)
  # require(future)
  # options(future.globals.maxSize = 30000 * 1024^2)
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
  # return(seurat)
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
  
  # require(future)
  # options(future.globals.maxSize = 65000 * 1024^2)
  # plan("multiprocess", workers = 2)
  
  anchors <- FindIntegrationAnchors(object.list = samples, anchor.features = features,
                                    normalization.method = "SCT", dims = 1:nCCA)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  print('done with integrate data')
  return(integrated)
}





analyze.seurat <- function(seurat, pcs.to.compute=200, pcs.to.analyze=40, debatch=FALSE) {
  #this function can adjust for different sample size
  require(Seurat)
  
  if(dim(seurat)[2] < 50){
    pcs.to.analyze=20
    pcs.to.compute=20
  }
  
  #skip this step if already DEBATCH, since SCT already run in debatching steps
  if(debatch==FALSE){
    DefaultAssay(seurat) <- "RNA"
    seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt", "percent.rp"),conserve.memory = T)
  }
  
  
  seurat <- RunPCA(seurat, npcs=pcs.to.compute)
  seurat <- FindNeighbors(seurat, dims=1:pcs.to.analyze)
  seurat <- FindClusters(seurat)
  
  return(seurat)}




seurat.SCT <- function(seurat, debatch=FALSE) {
  #this function can adjust for different sample size
  require(Seurat)
  
  #skip this step if already DEBATCH, since SCT already run in debatching steps
  if(debatch==FALSE){
    DefaultAssay(seurat) <- "RNA"
    seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt", "percent.rp"),conserve.memory = T)
  }
  
  if(debatch==FALSE){
    #if number of cells large 
    if(dim(seurat)[2] > 500000){
      print(paste0("Since number of cells large: ", 
                   dim(seurat)[2], " will use 'return.only.var.genes = TRUE' in SCTransform"))
      DefaultAssay(seurat) <- "RNA"
      # seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt", "percent.rp"),
      #                         conserve.memory = TRUE, 
      #                         return.only.var.genes = TRUE, 
      #                         verbose = TRUE)
      seurat <- SCTransform(seurat, vars.to.regress = c("percent.mt", "percent.rp"),
                            method = "glmGamPoi",
                            conserve.memory = TRUE, 
                            return.only.var.genes = TRUE, 
                            verbose = TRUE)
      
      
    }
  }
  return(seurat)
}





seurat.PCA <- function(seurat, pcs.to.compute = 200){
  if(dim(seurat)[2] < 50){
    pcs.to.compute=20
  }
  print(paste0("npcs to compute: ", pcs.to.compute))
  seurat <- RunPCA(seurat, npcs=pcs.to.compute)
  return(seurat)
}

seurat.findNeighbors <- function(seurat, pcs.to.analyze=40){
  if(dim(seurat)[2] < 50){
    pcs.to.analyze=20
  }
  seurat <- FindNeighbors(seurat, dims=1:pcs.to.analyze)
  return(seurat)
}

seurat.findClusters <- function(seurat){
  seurat <- FindClusters(seurat)
  return(seurat)
}

seurat.visualization <- function(seurat, pcs.to.analyze = 40){
  if(dim(seurat)[2] < 50){
    pcs.to.analyze=20
  }
  
  seurat <- RunUMAP(seurat, dims = 1:pcs.to.analyze, reduction.name = 'umap.rna', 
                    reduction.key = 'rnaUMAP_')
  
  if(dim(seurat)[2] <= 100){
    seurat <- RunTSNE(seurat, dims = 1:pcs.to.analyze, reduction.name="tsne.rna", 
                      reduction.key = "rnaTSNE_", perplexity=10)
  }
  if(dim(seurat)[2] > 100){
    seurat <- RunTSNE(seurat, dims = 1:pcs.to.analyze, reduction.name="tsne.rna", 
                      reduction.key = "rnaTSNE_", check_duplicates = FALSE)
    #duplicates here will be cells with the same PCA scores for the specified dimensions.
  }
  return(seurat)
}


analyze.seurat_2nd.stage <- function(seurat, pcs.to.compute=200, pcs.to.analyze=40) {
  #this function can adjust for different sample size
  require(Seurat)
  require(dplyr)
  
  if( (dim(seurat)[2] >=100) & (dim(seurat)[2] <=200)){
    pcs.to.analyze=20
    pcs.to.compute=20
  }
  #npcs = min(ncol(seurat), nrow(seurat)) - 1
  print(paste0("npcs to compute: ", pcs.to.compute))
  seurat <- seurat %>% RunPCA(npcs = pcs.to.compute)
  # approx=TRUE)  #use truncated SVD for faster computing
  seurat <- seurat %>% FindNeighbors(dims=1:pcs.to.analyze)
  seurat <- seurat %>% FindClusters()
  seurat <- seurat %>% RunUMAP(dims = 1:pcs.to.analyze, reduction.name = 'umap.rna', 
                               reduction.key = 'rnaUMAP_')
  
  #wont apply, since already filter seurat with with >=100 cells
  # if(dim(seurat)[2] <= 100){
  #   seurat <- seurat %>% RunTSNE(dims = 1:pcs.to.analyze, reduction.name="tsne.rna", 
  #                                reduction.key = "rnaTSNE_", perplexity=10)
  #   #higher perplexity, more smaller clusters
  #   #low perplexity, less clusters, and these clusters are bigger and closer together
  # }
  # if(dim(seurat)[2] > 100){
  seurat <- seurat %>% RunTSNE(dims = 1:pcs.to.analyze, reduction.name="tsne.rna", 
                               reduction.key = "rnaTSNE_", check_duplicates = FALSE)
  #duplicates here will be cells with the same PCA scores for the specified dimensions.
  
  
  return(seurat)
}

do.clusterings <- function(seurat, assay) {
  DefaultAssay(seurat) <- assay
  
  if(dim(seurat)[2] > 5000){
    resolutions <-  c(0.01, 0.05, 0.08, 0.1, 0.2)
  }
  if(dim(seurat)[2] <=5000){
    resolutions <- c(0.01, 0.02, 0.03, 0.05, 0.08, 0.1)
  }
  
  for(resolution in resolutions) {
    seurat <- FindClusters(seurat, resolution = resolution)
    seurat[[paste0("clusters_",assay, "_",resolution)]] = seurat$seurat_clusters
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
  #reference panel names:"dice", "monaco", "novershtern","blueprint","hpca","immgen", "mousernaseq"
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

get.covs <- function(seurat, preds, species, covs.clinical, curated.cts) {
  ref.by.species=list(Mouse=c("immgen", "mousernaseq"), Human=c("blueprint","hpca", "dice", "monaco","novershtern"))
  #this function is to get 10 cols of cell type prediction from singleR
  preds.labels <- preds.to.labels(
    subset.preds(preds,colnames(seurat))[,ref.by.species[[species]]])
  
  
  covs <- cbind(covs.clinical, preds.labels, default.cts = curated.cts)
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

get.covs_2ndstage<- function(covs, layer, curated.cell.types){
  covs_layer <- covs[covs$Celltype..malignancy. %in% layer,]
  covs_layer[[paste0("default.cts.", layer)]] <- curated.cell.types
  return(covs_layer)
}

get.Tisch.covs.clinical <- function(data.path, dataset, seurat.clusterings){ 
  data.path <- paste0(data.path, '/', dataset)
  
  #create meta.data file
  meta.file <- list.files(data.path, pattern =".*_CellMetainfo_table.tsv")
  meta.data <- read.table(file = paste0(data.path, '/', meta.file), 
                          sep = '\t', header = TRUE)
  #remove pre-computed umap from Tisch
  meta.data <- meta.data[,!colnames(meta.data) %in% c("UMAP_1","UMAP_2")]
  #filter cells to match with seurat after qc
  meta.data <- meta.data[meta.data$Cell %in% colnames(seurat.clusterings), ]
  return(meta.data)
}

###############

find.markers.by.clustering <- function(seurat, clusters) {
  DefaultAssay(seurat) <- "SCT"
  Idents(seurat) <- clusters
  
  #enable parallelization
  # require(future)
  # options(future.globals.maxSize = 64000 * 1024^2)
  # plan("multisession", workers = 2)
  
  markers <- FindAllMarkers(seurat)
  markers$diff <- markers$pct.1 - markers$pct.2
  markers <- markers %>% arrange(desc(diff))
  
  write.csv(markers, file = paste0("markers_", clusters, ".csv"))
  return(markers)
}


find.markers.by.all.clusterings <- function(seurat) {
  clusterings <- grep("clusters_.*", colnames(seurat@meta.data), value=TRUE)
  #check if clustering only has 1 group
  ret <- lapply(clusterings, function(clustering) {
    if (length(names(table(seurat[[clustering]]))) >= 2){
      find.markers.by.clustering(seurat,clustering)
    }
    else if(length(names(table(seurat[[clustering]]))) ==1){
      print(paste0("this clustering:  ", clustering, " only has 1 group"))
    }
  })
  names(ret) <- clusterings
  return(ret)
}



#####################################
###write all marker to xlsx 

write_all.markers <- function(fileName, all.markers){
  library(openxlsx)
  #check if all.markers file not empty
  if(class(all.markers[[1]])=="data.frame"){
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
  
}
####################

find.markers.Tisch.major.cts <- function(seurat, project){
  clusterings <- grep("Celltype.*", colnames(seurat@meta.data), value=TRUE)
  ret <- lapply(clusterings, function(clustering) {
    find.markers.by.clustering(seurat,clustering)
  })
  names(ret) <- clusterings
  write_all.markers(paste0(project, ".all.markers.Tisch.major.cts.xlsx"), ret)
  return(ret) 
}

#########################
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
  clusterings <- grep("clusters_.*", colnames(seurat@meta.data), value=TRUE)
  ret <- lapply(clusterings, function(clustering) {
    get.cell.and.cluster.stats.for.clustering(seurat, preds,clustering)
  })
  names(ret) <- clusterings
  return(ret)
}


write_cells.clusters.stats <- function(fileName, cell.and.cluster.stats){
  require(openxlsx)
  excel <- createWorkbook(fileName)
  
  #write multiple tabs/worksheets
  clusterings <- names(cell.and.cluster.stats)
  
  sapply(clusterings, function(clustering) {
    cell.stats <- cell.and.cluster.stats[[clustering]][["cell.stats"]]
    cluster.stats <- cell.and.cluster.stats[[clustering]][["cluster.stats"]]
    cluster_name <- gsub("clusters_", "", clustering)
    
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
  clustering = names(cell.and.cluster.stats)[1]
  cell.stats <- cell.and.cluster.stats[[clustering]][["cell.stats"]]
  cluster.stats <- cell.and.cluster.stats[[clustering]][["cluster.stats"]]
  
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


###########################

heatmap_artifacts_from_seurat_V5 <- function(seurat, clustering, markers, assay, min_lfg=0.0,
                                             max_p_adj = 0.05, limit=10, mask_mt_rp=TRUE) {
  Seurat::DefaultAssay(seurat) <- assay
  if(is.null(seurat[[assay]]@data)) {
    seurat <- Seurat::NormalizeData(seurat)
  }
  if(is.null(seurat[[assay]]@scale.data)) {
    seurat <- Seurat::ScaleData(seurat)
  }
  Seurat::Idents(seurat) <- clustering
  
  cluster_averages <- Seurat::AverageExpression(seurat, return.seurat = TRUE, assays = assay)
  
  
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
  
  extract_data_frame <- function(data_source, features, cluster_averages) {
    library(tibble)
    colnames(data_source) <- Idents(cluster_averages)
    data_source <- as.data.frame(data_source)
    data_source <- cbind(feature.name = features, data_source)
    rownames(data_source) <- data_source$feature.name
    return(data_source)
  }
  
  all_markers_lognorm <- extract_data_frame(cluster_averages[[cluster_averages@active.assay]]@layers[['data']], rownames(cluster_averages), cluster_averages)
  all_markers_scaled <- extract_data_frame(cluster_averages[[cluster_averages@active.assay]]@layers[['scale.data']], rownames(cluster_averages), cluster_averages)
  heatmap_markers_lognorm <- all_markers_lognorm[match(markers_use, all_markers_lognorm$feature.name), ]
  heatmap_markers_scaled <- all_markers_scaled[match(markers_use, all_markers_scaled$feature.name),]
  
  list(
    markers = markers %>% dplyr::mutate(cluster = as.character(.data$cluster)),
    all_markers_lognorm = all_markers_lognorm,
    all_markers_scaled = all_markers_scaled,
    heatmap_markers_lognorm = heatmap_markers_lognorm,
    heatmap_markers_scaled = heatmap_markers_scaled
  )
}

create.artifacts.ALL.clusterings <- function(seurat, all.markers){
  ##create artifacts for ALL clusterings
  clusterings_artifacts <- lapply(names(all.markers), function(clustering){
    print(clustering)
    #check if Idents has more than 1 level
    if (length(names(table(seurat[[clustering]]))) >= 2){
      artifacts <- heatmap_artifacts_from_seurat_V5(seurat=seurat, 
                                                    clustering = clustering, 
                                                    markers = all.markers[[clustering]], 
                                                    assay="SCT")
    }
  })
  names(clusterings_artifacts) <- names(all.markers)
  #remove empty list
  clusterings_artifacts <- clusterings_artifacts[unlist(lapply(clusterings_artifacts, length) != 0)]
  return(clusterings_artifacts)
}


prep_covs_Tisch <- function(covs){
  #this function is to (1) strip all non-clinical covs and only keep meta.data from cohort
  #(2) prep for input of artifacts in covs to write h5
  
  #remove cell barcodes in first column
  covs <- covs[,-1]
  
  qc_features <- c("nCount_RNA","nFeature_RNA", "percent.mt", "percent.rp", "nCount_SCT", "nFeature_SCT")
  
  clustering_names = grep("(Cluster)|(Celltype.*)",colnames(covs), value=T)
  
  
  covs.clinical <- covs[, !colnames(covs) %in% c(qc_features, clustering_names)]
  
  return(list(qc_features = qc_features, 
              clustering_names = clustering_names,
              extra_discrete_covs = covs.clinical))
}

update.seurat.list <- function(seurat.list, all.markers, tisch.markers){
  #special treatment for "all" layer
  #combine both tisch markers and all.markers from pipeline
  all.markers.and.tisch <- c(all.markers, tisch.markers)
  seurat.list$all$markers <- all.markers.and.tisch
  return(seurat.list)
}

create.artifacts.for.all.layers <- function(seurat.list.all.layers){
  #create artifacts for all layers
  layers <- lapply(names(seurat.list.all.layers), function(layer){
    print(paste0("working on layer: ", layer))
    seurat_layer <- seurat.list.all.layers[[layer]]$seurat
    covs_layer <- seurat.list.all.layers[[layer]]$covs
    markers_layer <- seurat.list.all.layers[[layer]]$markers
    
    #create artifact (covs) and clusterings for each layer:
    
    #prep for clusterings heatmap_artifacts
    clusterings_artifacts <- create.artifacts.ALL.clusterings(seurat_layer, markers_layer)
    
    #prep for covs
    covs.Tisch_layer <- prep_covs_Tisch(covs_layer)
    
    clustering_names.seurat <- grep("clusters_.*", colnames(seurat_layer@meta.data), value = TRUE)
    clustering_names <- c(clustering_names.seurat, covs.Tisch_layer$clustering_names)
    
    covs_artifacts_layer <- layer_artifacts_from_seurat(seurat_layer, 
                                                        qc_features=covs.Tisch_layer$qc_features, 
                                                        umap = "umap.rna", tsne = "tsne.rna",
                                                        clustering_names = clustering_names, 
                                                        extra_discrete_covs = covs.Tisch_layer$extra_discrete_covs,
                                                        extra_continuous_covs = NULL)
    
    
    #adding clustering artifacts into covs
    covs_artifacts_layer[["clusterings"]] <- clusterings_artifacts
    return(covs_artifacts_layer)
  })
  names(layers) <- names(seurat.list.all.layers)
  return(layers)
}

writing.iscvam.h5 <- function(seurat.list.all.layers, layers, dataset, fn){
  seurat <- seurat.list.all.layers$all$seurat
  # Call updated write_h5 that uses layer argument for Seurat v5
  write_h5(fn, seurat, layers, assays = c("RNA"))
}






#########################
##############

library(dplyr)
do.main.categories.lymphoids.pipeline <- function(cts){
  case_when(
    grepl("^NK.*", cts) ~ "NK cells",
    grepl("^B$", cts) ~ "B cells",
    
    grepl("^T.*", cts) ~ "T-cells",
    grepl("^CD4", cts) ~ "CD4",
    grepl("^CD8", cts) ~ "CD8",
    grepl("pDC", cts) ~ "pDC", 
    TRUE ~ cts)}

do.lymphoids <- function(cts){
  case_when(
    grepl("NK cells", cts) ~ "Lymphoids", 
    grepl("B cells", cts) ~ "Lymphoids", 
    grepl("T-cells", cts) ~ "Lymphoids", 
    grepl("CD4", cts) ~ "Lymphoids", 
    grepl("CD8", cts) ~ "Lymphoids", 
    grepl("pDC", cts) ~ "Lymphoids", 
    TRUE ~ cts)}

do.main.categories.myloids.pipeline <- function(cts){
  case_when(
    grepl("Mono.*", cts) ~ "Myloids", 
    grepl("Microglia.*", cts) ~ "Myloids", 
    TRUE ~ cts)}


do.maglinant <- function(cts){
  case_when(
    grepl(".*malignant.*", cts, ignore.case = TRUE) ~ "Malignant", 
    TRUE ~ cts)
} 

do.fibroblast <- function(cts){
  case_when(
    grepl("Fibroblast.*", cts) ~ "Fibroblasts", 
    TRUE ~ cts)
}

do.main.categories.pipeline <- function(cts){
  inverse <- c("Lymphoids", "Myloids", "Malignant", "Fibroblasts") 
  cts <- do.main.categories.lymphoids.pipeline(cts)
  cts <- do.lymphoids(cts)
  cts <- do.main.categories.myloids.pipeline(cts)
  cts <- do.maglinant(cts)
  cts <- do.fibroblast(cts)
  cts[!cts %in% inverse] <- "Others"
  
  return(cts)
}

###################

second_stage_analyze <- function(seurat, covs, layer){
  #this function perform all major steps for analyzing 2nd stage for 1 layer
  
  #subset cells
  subset.cts <- subset(seurat, idents = layer)
  print(paste0("number of cells of Seurat for this layer: ", dim(subset.cts)[2]))
  
  cols.to.remove <- grep("SCT_.*", colnames(subset.cts@meta.data), value=TRUE)
  subset.cts[[cols.to.remove]] <- NULL
  
  #perform clustering stage
  subset.seurat <- analyze.seurat_2nd.stage(subset.cts)
  subset.clusterings <- do.clusterings(subset.seurat, "SCT")
  
  layer_all.markers <- find.markers.by.all.clusterings(subset.clusterings)
  layer_export.rna.markers <- write_all.markers(paste0(layer, "_all.markers.xlsx"), layer_all.markers)
  
  layer_cell.and.cluster.stats <- get.cell.and.cluster.stats.for.all.clusterings(subset.clusterings, singler)
  export.cells.clusters.stats <- write_cells.clusters.stats(paste0(layer, '_cells.and.clusterings.xlsx'), layer_cell.and.cluster.stats)
  layer_curated.cell.types <- do.curated.cell.types(layer_cell.and.cluster.stats)                                                          
  covs_layer <- get.covs_2ndstage(covs, layer, layer_curated.cell.types)
  
  return(list(seurat = subset.clusterings, 
              covs = covs_layer, 
              markers = layer_all.markers))
}


do.multiStages <- function(seurat.list){
  #prep for 2nd stage
  #extract covs and seurat
  covs <- seurat.list$all$covs
  seurat <- seurat.list$all$seurat
  
  #renaming main categories based on major cell types in Tisch
  main.categories <- covs$Celltype..malignancy.
  print("number of cells for each layer: ")
  print(table(covs$Celltype..malignancy.))
  
  #update main.categories to seurat and covs
  seurat$main.categories <- main.categories
  Idents(seurat) <- seurat$main.categories
  
  #extract names of main.categories
  dataset.main.categories <- names(table(covs$Celltype..malignancy.))
  
  ####
  #check if number of cells in each layer >= 100 cells
  check.numb.cells <- janitor::tabyl(covs, Celltype..malignancy.) %>% 
    mutate(do.secondStage = ifelse(n >=100, TRUE, FALSE))
  print(check.numb.cells)
  # 
  dataset.main.categories <- check.numb.cells[check.numb.cells$do.secondStage==TRUE, ]$Celltype..malignancy.
  
  print(paste0("names of main.categories: ", dataset.main.categories))
  
  #perform 2nd stage for all layers
  all.second.stages <- lapply(dataset.main.categories, function(layer){
    print(paste0("working on layer: ", layer))
    second.stage.lst <- second_stage_analyze(seurat, covs, layer)
    save(second.stage.lst, file = paste0(layer, "_seurat.list.RData"))
    return(second.stage.lst)
  })
  names(all.second.stages) <- dataset.main.categories
  
  #combine 1st stage and 2nd stage
  seurat.list.all.layers <- c(seurat.list, all.second.stages)
  return(seurat.list.all.layers)
  
}

