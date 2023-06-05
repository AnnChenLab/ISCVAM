library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)

library(BSgenome.Hsapiens.UCSC.hg38)
library(rhdf5)
library(dplyr)
library(SingleR)
library(future)


set.seed(1234)

create_seurat_multiome <- function(data.path, meta.data=FALSE){
  #this Seurat is without meta.data
  
  if(any(endsWith(list.files(data.path), ".h5"))){
    path.10x <- paste0(data.path, "/", list.files(data.path, pattern = ".*.h5"))
    print(path.10x)
    counts = Read10X_h5(path.10x)
  }else{
    path.10x <- paste0(data.path, "/filtered_feature_bc_matrix")
    print(path.10x)
    counts <- Read10X(path.10x)
  }
  
  fragpath <- list.files(data.path, pattern='.*_fragments.tsv.gz', recursive = F)
  
  library(EnsDb.Hsapiens.v86)
  # get gene annotations for hg38
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  
  # create a Seurat object containing the RNA adata
  seurat <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )
  
  # create ATAC assay and add it to the object
  seurat[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0(data.path,'/', fragpath[1]),
    annotation = annotation
  )
  if (meta.data==TRUE){
    seurat <- create_seurat.w.meta.data(data.path, seurat)
    return(seurat)
  }else{
    return(seurat)
  }
}

#############
#function to create seurat.multiome for many samples

create.seurat.raw.batch <- function(sample.table, meta.data=meta.data){
  seurat.raw.batch <- lapply(sample.table$data.path, function(data.path){
    print(data.path)
    seurat.raw <- create_seurat_multiome(data.path, meta.data=meta.data)
  })
  names(seurat.raw.batch) <- sample.table$dataset.names
  return(seurat.raw.batch)
}

################################
# PEAK CALLING


create_peaks_assay <- function(peaks, seurat){
  library(GenomicRanges)
  DefaultAssay(seurat) <- "ATAC"
  #input 'peaks' can be either from 'peak calling' from MACS2 or from 'peaks.bed' file
  peaks <- makeGRangesFromDataFrame(df = peaks)
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions

  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  library(future)
  plan("multisession", workers = 10)
  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(seurat),
    features = peaks,
    cells = colnames(seurat)
  )
  
  seurat[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    #fragments = paste0(data.path,'/', fragpath[1])
    fragments = Fragments(seurat)
    annotation = Annotation(seurat)
  )
  return(seurat)
}

#only apply for LIST of seurat objects
get.peak.calling.batch<- function(sample.table, seurat.raw.ls){
  library(BiocParallel)
  
  data.path.ls <- sample.table$data.path
  ### implement peaks calling for all seurat objects
  peaks.calling_ls <- lapply(data.path.ls, function(path){
    #extract dataset name
   data.name <- tail(strsplit(path, "/")[[1]], 1)
    
   #function"PEAKS_CALLING" will vary based on projects, thus find this function in prep_xxx.R
   #some projects already provided BED file like ovarian, some not, thus need to run MACS2 
    seurat <- peaks_calling(path, seurat.raw.ls[[data.name]])
    return(seurat)
  })
  names(peaks.calling_ls) <- sample.table$dataset.names
  return(peaks.calling_ls)
}


#############
### check dimensions of individual datasets before and after QC ATAC and RNA

retrieve_dimension <- function(seurat){
  dimensions <- sapply(names(seurat@assays), function(assay){
    dim(seurat[[assay]]) })
  df <- data.frame(dimensions)
  df_1 <- cbind(df[1,], nCells = df[2, 1])
  return(df_1)
}


retrieve_dimension_seurat_batch <- function(list_seurat){
  dimensions <- lapply(list_seurat, retrieve_dimension)
  dimensions <- as.data.frame(do.call(rbind, lapply(dimensions, as.data.frame)))
  write.csv(dimensions, file = "track_modality_dimensions_before_qc_batch.csv")
  return(dimensions)
}


################
#UNION of all peaks in all datasets

create.union.seurat <- function(seurat_raw_batch_peaks){
  #create a UNION/combining all peaks of 8 samples (comparing "intersect" peaks)'''
  library(BiocParallel)
  #seurat_raw_batch_peaks = seurat.raw.ls 
  #ADD DATASET'S NAME TO 8 SAMPLES, PREPARE TO MERGE
  seurat_batch_update <- bplapply(names(seurat_raw_batch_peaks), function(sample){
    seurat_raw_batch_peaks[[sample]]@meta.data$dataset <- sample
    return(seurat_raw_batch_peaks[[sample]])
  })
  names(seurat_batch_update) <- names(seurat_raw_batch_peaks)


  #MERGE (UNION) 8 SAMPLES (if default, 'all' = FALSE meaning intersect)
  union_seurat_8samples <- merge(
    x = seurat_batch_update[[1]],
    y = seurat_batch_update[-1],
   # add.cell.ids = names(seurat_batch_update),
    all=TRUE #UNION
    )

   return(union_seurat_8samples)
}

###########
do.singler <- function(raw.reads, single.ref.rds) {
  require(SingleR)
  require(future)
  load(single.ref.rds)
  options(future.globals.maxSize = 800000 * 1024^2)
  preds <- sapply(names(refs), function(ref) {
    sapply(c("main","fine"), function(r) {
      SingleR(raw.reads, refs[[ref]], refs[[ref]][[paste0("label.",r)]] ) }
    )
  })
}



#############
#QC in ATAC
qc_ATAC <- function(seurat, meta.data =TRUE){
  options(future.globals.maxSize = 800000 * 1024^2)
  DefaultAssay(seurat) <- "ATAC"
  
  seurat <- NucleosomeSignal(seurat)
  seurat <- TSSEnrichment(seurat)
  
  if (meta.data==TRUE){
    seurat$pct_reads_in_peaks <- seurat$atac_peak_region_fragments / seurat$atac_fragments * 100
    
    seurat <- subset(
      x = seurat,
      subset = atac_peak_region_fragments > 3000 &
        atac_peak_region_fragments < 20000 &
        pct_reads_in_peaks > 15 &
        nucleosome_signal < 4 &
        TSS.enrichment > 2
    )
  }else{ #if meta.data not avai, only QC ATAC based on nCount_ATAC
    seurat <- subset(
      x = seurat,
      subset = nCount_ATAC < 100000 &
        nCount_RNA < 25000 &
        nCount_ATAC > 1000 &
        nCount_RNA > 1000 &
        nucleosome_signal < 2 &
        TSS.enrichment > 1
    )
  }
  return(seurat)
}

###############
## QC in RNA
qc_RNA <- function(seurat, species) {
  DefaultAssay(seurat) <- "RNA"
  require(dplyr)
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
  seurat[,seurat[["percent.mt"]] <= 20 & seurat[["nFeature_RNA"]] >= 500]
  
}

####################
qc_rna_atac_batch <- function(seurat_raw_batch_peaks, meta.data=meta.data){
  # QC individual seurat objects before merging 
  seurat_qc_atac_rna_batch <- lapply(seurat_raw_batch_peaks, function(seurat){
    seurat.qc.atac <- qc_ATAC(seurat, meta.data = meta.data) 
    seurat.qc.rna <- qc_RNA(seurat.qc.atac, "Human")
    return(seurat.qc.rna)
  })
  names(seurat_qc_atac_rna_batch) <- names(seurat_raw_batch_peaks)
  return(seurat_qc_atac_rna_batch)
}



#################
integrate.seurat.by.RNA <- function(seurat.merge, nGenes, nCCA, split.by="patient"){ 
    # nCCA=30
  # nGenes=6000
  library(BiocParallel)
  
  patients <- SplitObject(seurat.merge, split.by = split.by) 
  
 
  #change default assay in each individual seurat
  patients <- lapply(patients, function(seurat){
                    DefaultAssay(seurat) <- "RNA"
                    return(seurat)}) %>% 
    lapply(SCTransform, vars.to.regress = "percent.mt", 
                        method = "glmGamPoi")
    #here SCTransform use "glmGamPoi" instead of "poisson" because pt "012320M" not working w "poisson"
  #also "glmGamPoi" was recommended in this paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9
  
  features <- SelectIntegrationFeatures(object.list = patients, 
                                        anchorfvf.nfeatures = 6000, # nfeatures for FindVariableFeatures. Used if VariableFeatures have not been set for any object in object.list.
                                        nfeatures = nGenes)
  
  patients <- PrepSCTIntegration(object.list = patients, anchor.features = features)
  
  require(future)
  options(future.globals.maxSize = 800000 * 1024^2)
  plan("multisession", workers = 6)
  
  anchors <- FindIntegrationAnchors(object.list = patients, anchor.features = features,
                                    normalization.method = "SCT", dims = 1:nCCA)
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  print('done with integrate data')
  return(integrated)
}

##########
### analyze RNA

analyze_rna <- function(seurat,  pcs.to.compute=200, pcs.to.analyze=40, debatch=FALSE) {
  
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
  seurat <- seurat %>% RunUMAP(dims = 1:pcs.to.analyze,  reduction.name = 'umap.rna', 
                               reduction.key = 'rnaUMAP_')
  seurat <- seurat %>% RunTSNE(dims = 1:pcs.to.analyze, reduction.name="tsne.rna", 
                               reduction.key = "rnaTSNE_")
}


do.clusterings <- function(seurat, assay) {
  DefaultAssay(seurat) <- assay
  for(resolution in c(0.6,0.8, 1, 1.2, 2, 4)) {
    seurat <- FindClusters(seurat, resolution = resolution)
    seurat[[paste0("clusters_",assay, "_",resolution)]] = seurat$seurat_clusters
  }
  return(seurat)
}


do.clusterings.WNN <- function(seurat) {
  #no need to specify assay since WNN is already combining rna and atac
  #DefaultAssay(seurat) <- assay
  
  for(resolution in c(0.6,0.8, 1, 1.2, 2, 4)) {
    seurat <- FindClusters(seurat, graph.name = "wsnn", algorithm = 3, verbose = FALSE, 
                 resolution=resolution)
    seurat[[paste0("clusters_wsnn_",resolution)]] = seurat$seurat_clusters
  }
  return(seurat)
}

#############
#analyze ATAC


analyze_atac <- function(RWPE.combined, debatch=FALSE){

  DefaultAssay(RWPE.combined) <- "peaks"

  RWPE.combined <- RunTFIDF(RWPE.combined)
  RWPE.combined <- FindTopFeatures(RWPE.combined, min.cutoff = 20) #Only peaks presented in more than 20 cells were selected
  RWPE.combined <- RunSVD(RWPE.combined)
  
  if(debatch){
    RWPE.combined <- RunUMAP(RWPE.combined, reduction = 'harmony', dims = 2:50, reduction.name="umap.harmony", reduction.key="harmonyUMAP_")
    RWPE.combined <- RunTSNE(RWPE.combined, reduction = "harmony", dims = 2:50, reduction.name="tsne.harmony", reduction.key="harmonyTSNE_")
    RWPE.combined <- FindNeighbors(RWPE.combined, reduction = 'harmony', dims = 2:50)
  } 
  else{
    RWPE.combined <- RunUMAP(RWPE.combined, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    RWPE.combined <- RunTSNE(RWPE.combined, reduction = 'lsi', dims = 2:50, reduction.name = "tsne.atac", reduction.key = "atacTSNE_")
    RWPE.combined <- FindNeighbors(RWPE.combined, reduction = 'lsi', dims = 2:50)
    }
  return(RWPE.combined)
  
}

###################
# debatch ATAC using HARMONY

integrate.seurat.by.ATAC <- function(union_peaks){
  #https://www.10xgenomics.com/resources/analysis-guides/batch-effect-correction-in-chromium-single-cell-atac-data
  library(harmony)
  
  merge_atac_analyzed <- analyze_atac(union_peaks, debatch=FALSE)
  
  hm.integrated <- RunHarmony(object = merge_atac_analyzed, group.by.vars = 'dataset', 
                              reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
  
  hm.integrated <- analyze_atac(hm.integrated, debatch=TRUE)
  return(hm.integrated)
  
}


###############
### WNN
# ###########################
# #Joint UMAP visualization
# # build a joint neighbor graph using both assays

build_WNN <- function(seurat.analyzed.atac, debatch.atac=FALSE){
  if (debatch.atac ==TRUE){
    reduction.list = list("pca", "harmony")
  }
  else{ #meaning debatch.atac=FALSE
    reduction.list = list("pca", "lsi")
  }
  
  seurat.analyzed.atac <- FindMultiModalNeighbors(
  object = seurat.analyzed.atac,
  reduction.list = reduction.list, 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE)

  # # build a joint UMAP visualization
  seurat.analyzed.atac <- RunUMAP(
    object = seurat.analyzed.atac,
    nn.name = "weighted.nn",
    verbose = TRUE, 
    reduction.name = "wnn.umap", 
    reduction.key = "wnnUMAP_"
  )

  seurat.analyzed.atac <- RunTSNE(
    object = seurat.analyzed.atac,
    nn.name = "weighted.nn",
    verbose = TRUE, 
    reduction.name = "wnn.tsne", 
    reduction.key = "wnnTSNE_"
  )
  return(seurat.analyzed.atac)
}


###############
#####################
###### TF #############

TF_analysis <- function(RWPE.combined.analyzed){
  library(motifmatchr)
  library(JASPAR2020)
  library(TFBSTools)
  
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(EnsDb.Hsapiens.v86)
  
  # Get a list of motif position frequency matrices from the JASPAR database
    pfm <- getMatrixSet(
      x = JASPAR2020,
      opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE, 
                  species = 'Homo sapiens')
    )

    DefaultAssay(RWPE.combined.analyzed) <- "peaks"
    # add motif information
    RWPE.combined.analyzed <- AddMotifs(
      object = RWPE.combined.analyzed,
      genome = BSgenome.Hsapiens.UCSC.hg38,
      pfm = pfm
    )


    RWPE.combined.analyzed <- RunChromVAR(
      object = RWPE.combined.analyzed,
      genome = BSgenome.Hsapiens.UCSC.hg38
    )
    return(RWPE.combined.analyzed)

}

###############
find.markers.by.clustering <- function(seurat, clusters, assay) {
  #REF: https://github.com/satijalab/seurat/issues/1717 
  if ((assay=="SCT") || (assay =="integrated")){
    #when doing clustering, clustering based on SCT, which store < 36k genes 
    #but when finding marker genes, must use RNA assay which store 36k genes
    
    #INPUT assay here is "SCT" because using SCT for clustering, thus clustering column name included "SCT" word
    assay <- "RNA"
    DefaultAssay(seurat) <- "RNA"
  }else{
    DefaultAssay(seurat) <- assay
  }
  
  Idents(seurat) <- clusters
  
  #enable parallelization
  require(future)
  options(future.globals.maxSize = 800000 * 1024^2)
  plan("multisession", workers = 4)
  
  print(paste0("assay of markers is: ", DefaultAssay(seurat)))
  markers <- FindAllMarkers(seurat, assay=assay)
  markers$diff <- markers$pct.1 - markers$pct.2
  return(markers)
}


find.markers.by.all.clustering <- function(seurat, group=1, assay) {
  if (group ==1){
    clusterings <- c(paste0("clusters_", assay, "_0.6"), 
                     paste0("clusters_", assay, "_0.8"),
                     paste0("clusters_", assay, "_1") 
                     )
  }
  if (group ==2){
    clusterings <- c(paste0("clusters_", assay, "_1.2"),
                      paste0("clusters_", assay, "_2"),
                     paste0("clusters_", assay, "_4"))
  }
  if (group ==0){
    clusterings <- c(paste0("clusters_", assay, "_0.6"),
                    paste0("clusters_", assay, "_0.8"),
                    paste0("clusters_", assay, "_1"),
                    paste0("clusters_", assay, "_1.2"), 
                    paste0("clusters_", assay, "_2"), 
                    paste0("clusters_", assay, "_4"))
  }
  if(group=="wnn"){
    clusterings <- grep("clusters_wsnn_.*", colnames(seurat@meta.data), value = TRUE)
  }

  ret <- lapply(clusterings, function(clustering) {
    print(paste0('start with ', clustering))
    find.markers.by.clustering(seurat,  seurat[[clustering]], assay)
    
  })
  names(ret) <- clusterings
  ret
}



##############################
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

###################

plot.heatmap.seurat <- function(seurat, markers, fileName, assay, marker.set.name) {
  library(Seurat)
  pdf(file= fileName, height=28,width=20)
  
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat, assay = "RNA")
  seurat <- ScaleData(seurat)
  
  DefaultAssay(seurat) <- "peaks"
  seurat <- NormalizeData(seurat, assay = "peaks")
  seurat <- ScaleData(seurat)
  
  
  DefaultAssay(seurat) <- assay
  cluster.averages <- AverageExpression(seurat, return.seurat = TRUE)
  
  DefaultAssay(cluster.averages) <- assay
  cluster.averages <- NormalizeData(cluster.averages, assay = assay)
  cluster.averages <- ScaleData(cluster.averages)
  
  require(ggplot2)
  require(dplyr)
  min.lfg <- 0.0
  max.p.adj <- 0.05
  #limit <- 20
  n.clusters <- length(unique(markers$cluster))
 limit <- 140/n.clusters
  
  for(mask.mt.rp in c(T,F)) {
    markers.use=subset(markers,avg_log2FC > min.lfg & 
                         p_val_adj < max.p.adj & 
                         !(mask.mt.rp & (grepl("^RP[SL]", gene) | grepl("^MT", gene)))) %>%
      group_by(cluster) %>% 
      top_n(limit, -p_val_adj) %>% 
      top_n(limit, abs(pct.1-pct.2)) %>%
      #if clusters is in FACTOR, then "arrange(as.factor(cluster))"
      #because as.character can't sort alpha_numeric (cant sort 10 instead of 1)
      arrange(as.factor(cluster))
    # arrange(as.character(cluster))
    
    clus.avg_list <- export.cluster.avg.and.heatmap(cluster.averages, markers.use,
                                   xlsx.name=paste0('clus.avg.', marker.set.name, '.xlsx'), assay)
    
    markers.use <- as.character(markers.use$gene)
    
    print(DoHeatmap(cluster.averages, features = markers.use,
                    draw.lines = FALSE) + 
            ggtitle(paste0("Average expression; mask.mt.rp=", mask.mt.rp)))
    
    print(DoHeatmap(seurat, features = markers.use, assay=assay) +
            ggtitle(paste0("test=wilcox; zscore; mask.mt.rp=", mask.mt.rp)))
    
    print(DoHeatmap(seurat, features = markers.use, slot="data") +
            ggtitle(paste0("test=wilcox; log normalized; mask.mt.rp=", mask.mt.rp)))
    
    if(nrow(seurat@assays[[seurat@active.assay]]@counts) > 0) {
      print(DoHeatmap(seurat, features = markers.use, slot="counts" ) +
              ggtitle(paste0("test=wilcox; counts; mask.mt.rp=", mask.mt.rp)))
    }
    
  }
  dev.off()
  return(clus.avg_list)
}

####################
plot.heatmap_ALL <- function(marker.sets, seurat, assay){
  clus.avg_ALL <- lapply(names(marker.sets), function(marker.set.name){
  print(marker.set.name)
  markers = marker.sets[[marker.set.name]]

  Idents(seurat) = seurat[[marker.set.name]]
  
  fileName <- paste0("report.heatmap.",marker.set.name, ".pdf") 
  
  clus.avg_list <- plot.heatmap.seurat(seurat, markers, fileName, assay = assay, marker.set.name)
  return(clus.avg_list)
})
  names(clus.avg_ALL) <- names(marker.sets)
  return(clus.avg_ALL)
  
}

#######################

###############################
do.curated.cell.types <- function(cell.and.cluster.stats, clustering) {
  #cell.stats <- cell.and.cluster.stats[["cell.stats", clustering]]
  #cluster.stats <- cell.and.cluster.stats[["cluster.stats", clustering]]
  
  cell.stats <- cell.and.cluster.stats[[clustering]]$cell.stats
  cluster.stats <- cell.and.cluster.stats[[clustering]]$cluster.stats
  
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

#############################
#read results of singleR 

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

get.covs <- function(seurat, preds, species, default.ct) {
  ref.by.species=list(Mouse=c("immgen", "mousernaseq"), Human=c("blueprint","hpca", "dice", "monaco","novershtern"))
  #this function is to get 10 cols of cell type prediction from singleR
  preds.labels <- preds.to.labels(subset.preds(preds,colnames(seurat))[,ref.by.species[[species]]])

  preds.labels <- data.frame(preds.labels)
  covs <- cbind(default.ct = default.ct, preds.labels)
  return(covs)
}


################


get.cell.and.cluster.stats.for.clustering <- function(seurat, preds, clustering="clusters_SCT_0.6") {
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
    group_by(clusterN) %>%
    mutate(freq=round(n/sum(n), digits = 2))
  
  majority.stats <- cluster.stats %>% group_by(clusterN) %>% filter(freq==max(freq))
  list(cell.stats=cell.stats, cluster.stats=majority.stats)
}

get.cell.and.cluster.stats.for.all.clusterings <- function(seurat, preds) {
  # clusterings <- c("seurat_clusters_0.6",
  #                  "seurat_clusters_0.8", 
  #                  "seurat_clusters_1", 
  #                  "seurat_clusters_1.2",
  #                  "seurat_clusters_2",
  #                  "seurat_clusters_4")
  # "clusters_infomap")
  
  clusterings <- colnames(seurat@meta.data)[grep("clusters_.*", colnames(seurat@meta.data))]

  res <- lapply(clusterings, function(clustering) {
    get.cell.and.cluster.stats.for.clustering(seurat, preds,clustering)
  })
  names(res) <- clusterings
  return(res)
  
}

write_cells.clusters.stats <- function(fileName, cell.and.cluster.stats, seurat){
  require(openxlsx)
  excel <- createWorkbook(fileName)
  
  clusterings <- colnames(seurat@meta.data)[grep("clusters_.*", colnames(seurat@meta.data))]
  
  
  sapply(clusterings, function(clustering) {
    cell.stats <- cell.and.cluster.stats[[clustering]][['cell.stats']]

    cluster.stats <- cell.and.cluster.stats[[clustering]][["cluster.stats"]]
    cluster_name <- gsub('clusters', '', clustering)
    
    # Add worksheets to workbook
    addWorksheet(excel, sheetName = paste0(cluster_name, '_cell.stats'))
    writeData(excel, sheet = paste0(cluster_name, '_cell.stats'), cell.stats, colNames = TRUE)
    
    addWorksheet(excel, sheetName = paste0(cluster_name, '_cluster.stats'))
    writeData(excel, sheet = paste0(cluster_name, '_cluster.stats'), cluster.stats, colNames=TRUE)
  })
  
  saveWorkbook(excel, file = fileName, overwrite = TRUE)
}


########################
#this is default.ct for rna version i.e. finding rna markers for default.ct
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

###############
#plot heatmap RNA for default.ct



plot.heatmap.default.ct_rna <- function(annotation.default.cell.types, seurat.wnn.clusterings){
  ###plot heatmap for rna markers for default.ct
  
  DefaultAssay(seurat.wnn.clusterings) <- "RNA"
  
  order.ct <- names(table(annotation.default.cell.types$markers$cluster))
  
  Idents(seurat.wnn.clusterings) <- seurat.wnn.clusterings$default.ct
  levels(seurat.wnn.clusterings) <- order.ct
  
  
  clus.avg_default.ct_rna <- plot.heatmap.seurat(seurat.wnn.clusterings, 
                                                 annotation.default.cell.types$markers, 
                                                 "heatmap_default.ct_RNA.markers.pdf", 
                                                 "RNA", "default.ct_RNA")
  return(clus.avg_default.ct_rna)
}


######################
#default.ct for peak version i.e. finding peak markers for default.ct + export clus.avg 
default.ct_find.peak.markers.and.plot.heatmap <- function(seurat.wnn.clusterings, default.ct,
                                                          annotation.default.cell.types){
  #find peak markers for default.ct
  seurat.wnn.clusterings$default.ct <- default.ct
  
  order.ct <- names(table(annotation.default.cell.types$markers$cluster))
  Idents(seurat.wnn.clusterings) <- seurat.wnn.clusterings$default.ct
  levels(seurat.wnn.clusterings) <- order.ct
  
  
  DefaultAssay(seurat.wnn.clusterings) <- "peaks"
  default.ct_peaks_markers <- find.markers.by.clustering(seurat.wnn.clusterings, "default.ct", assay="peaks")
  
  #reorder marker list w correct cluster order
  default.ct_peaks_markers_order <- default.ct_peaks_markers %>%  
    mutate(cluster = factor(default.ct_peaks_markers$cluster, 
                            levels =  order.ct)) %>% arrange(cluster)
  
  
  #plot heatmap for default.ct for peak assay
  clus.avg_default.ct_peaks <- plot.heatmap.seurat(seurat.wnn.clusterings,  
                                                   markers =default.ct_peaks_markers_order, 
                                                   fileName = "heatmap_default.ct_peak.markers.pdf", 
                                                   assay = "peaks", 
                                                   marker.set.name = "default.ct_peaks")
  
  return(list(clus.avg = clus.avg_default.ct_peaks, 
              markers =default.ct_peaks_markers_order ))
}


#############################
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
  print(DimPlot(seurat, reduction = "umap.rna"))
  print(DimPlot(seurat, reduction = "tsne.rna"))
  print(PCAPlot( seurat))
  
  dev.off()
  
  #######
  #prep for exporting cluster avg and heatmap for default.ct
  
  cluster.averages <- create.cluster.avg(seurat, curated.cell.types)
 
  require(dplyr)
  min.lfg <- 0.0
  max.p.adj <- 0.05
  limit <- 20
  mask.mt.rp <- T
  markers.use=subset(marker.sets$markers,avg_log2FC > min.lfg & 
                       p_val_adj < max.p.adj & 
                       !(mask.mt.rp & (grepl("^RP[SL]", gene) | grepl("^MT", gene)))) %>%
    group_by(cluster) %>% 
    top_n(limit, -p_val_adj) %>% 
    top_n(limit, abs(pct.1-pct.2)) %>%
    arrange(as.factor(cluster))
  
  xlsx.name <- "cluster.averages_default.ct.xlsx"
  export.cluster.avg.and.heatmap(cluster.averages, markers.use, xlsx.name)
  
  fileName <- paste0(version, ".report.heatmap.","default.ct", ".pdf")
  plot.heatmap.seurat(seurat, marker.sets$markers, fileName)
  
}
###################################
do.sub.categories <- function(cts){
    sub.categories <- dplyr::case_when(
      grepl("^B .*", cts) ~ "B cells", 
      
      grepl("^CD8+.*", cts) ~ "CD8",
      grepl("^CD4+.*", cts) ~ "CD4",
      
      grepl("^CD14 .*", cts) ~ "Myeloids",
      grepl("^CD16 .*", cts) ~ "Myeloids",
      grepl("ASDC", cts) ~ "Myeloids", 
      
      grepl("HSPC", cts) ~ "Others",
      grepl("Eryth", cts) ~ "Others",
      grepl("ILC", cts) ~ "Others",
      
      grepl("^NK.*", cts) ~ "NK cells",
      
      TRUE ~ cts
    )
    return(sub.categories)
    }

do.main.categories <- function(sub.categories){
  main.categories <- dplyr::case_when(
    grepl("^B .*", sub.categories) ~ "Lymphoids", 
    grepl("^NK.*", sub.categories) ~ "Lymphoids",
    
    grepl("gdT", sub.categories) ~ "Lymphoids",
    grepl("dnT", sub.categories) ~ "Lymphoids",
    
    grepl("CD8", sub.categories) ~ "Lymphoids",
    grepl("CD4", sub.categories) ~ "Lymphoids",
    grepl("Treg", sub.categories) ~ "Lymphoids",
    
    grepl("pDC", sub.categories) ~ "Myeloids", 
    grepl("(cDC1)|(cDC2)", sub.categories) ~ "Myeloids",
    
    grepl("Platelet", sub.categories) ~ "Others",
    grepl("Plasmablast", sub.categories) ~ "Others",
    
    TRUE ~ sub.categories)
  return(main.categories)
}

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


################
do.main.sub.categories <- function(default.ct, seurat.wnn){
  sub.categories <- do.sub.categories(default.ct)
  main.categories <- do.main.categories(sub.categories)
  
  seurat.wnn$default.ct <- default.ct
  seurat.wnn$sub.categories <- sub.categories
  seurat.wnn$main.categories <- main.categories
  
  rna.main.cat <- find.markers.by.clustering(seurat.wnn, "main.categories", "RNA")
  peaks.main.cat <- find.markers.by.clustering(seurat.wnn, "main.categories", "peaks")
  
  rna.sub.cat <- find.markers.by.clustering(seurat.wnn, "sub.categories", "RNA")
  peaks.sub.cat <- find.markers.by.clustering(seurat.wnn, "sub.categories", "peaks")
  
  require(openxlsx)
  excel <- createWorkbook("main.sub.categories.markers.xlsx")
  
  addWorksheet(excel, sheetName = "main.cat_rna.markers")
  writeData(excel, sheet = "main.cat_rna.markers", rna.main.cat, colNames=TRUE)
  
  addWorksheet(excel, sheetName = "main.cat_peak.markers")
  writeData(excel, sheet = "main.cat_peak.markers", peaks.main.cat, colNames=TRUE)
  
  addWorksheet(excel, sheetName = "sub.cat_rna.markers")
  writeData(excel, sheet = "sub.cat_rna.markers", rna.sub.cat, colNames=TRUE)
  
  addWorksheet(excel, sheetName = "sub.cat_peak.markers")
  writeData(excel, sheet = "sub.cat_peak.markers", peaks.sub.cat, colNames=TRUE)
  
  saveWorkbook(excel, file = "main.sub.categories.markers.xlsx", overwrite = TRUE)
  
  return(list(main.categories = list(main.categories = main.categories, 
                                     main.cat_rna.markers = rna.main.cat, 
                                     main.cat_peak.markers = peaks.main.cat
                                     ), 
              sub.categories = list(sub.categories = sub.categories, 
                                    sub.cat_rna.markers = rna.sub.cat, 
                                    sub.cat_peak.markers = peaks.sub.cat)))
}


################################

create.cluster.avg <- function(motif_analysis, clustering_resolution){
  Idents(motif_analysis) <- clustering_resolution
  DefaultAssay(motif_analysis) <- "RNA"
  
  motif_analysis <- NormalizeData(motif_analysis, assay = "RNA")
  motif_analysis <- ScaleData(motif_analysis)
  cluster.averages <- AverageExpression(motif_analysis, return.seurat = TRUE)
  
  
  DefaultAssay(cluster.averages) <- "SCT"
  cluster.averages <- NormalizeData(cluster.averages, assay = "SCT")
  cluster.averages <- ScaleData(cluster.averages)
  
  
  DefaultAssay(cluster.averages) <- "peaks"
  cluster.averages <- NormalizeData(cluster.averages, assay = "peaks")
  cluster.averages <- ScaleData(cluster.averages)
  
  
  return(cluster.averages)
}

# add.feature.names <- function(dataset){
#   library(tibble)
#   feature.names <- rownames(dataset)
#   dataset <- dataset %>% as_tibble %>% 
#     add_column(feature.names=feature.names, .before=1) %>% 
#     as.data.frame()
#  # rownames(dataset) <- dataset$feature.names #because rownames may have duplicates, BUT rownames of df not allow duplicates
#   return(dataset)
# }
  

export.cluster.avg.and.heatmap <- function(cluster.averages, markers.use, xlsx.name, assay){
  
  library(openxlsx)
  wb<-createWorkbook(xlsx.name)
  
  addWorksheet(wb, sheetName = "all_markers_lognorm")
  all_markers_lognorm <- as.matrix(cluster.averages[[assay]]@data) %>% add.feature.names
  writeData(wb, sheet = "all_markers_lognorm", all_markers_lognorm, rowNames=TRUE)
  
  
  addWorksheet(wb, sheetName = "heatmap_markers_lognorm")
  heatmap_markers_lognorm <- as.matrix(cluster.averages[[assay]]@data[markers.use$gene,]) %>% add.feature.names
  writeData(wb, sheet = "heatmap_markers_lognorm", heatmap_markers_lognorm, rowNames=TRUE)
  
  
  addWorksheet(wb, sheetName = "all_markers_scaled")
  all_markers_scaled <- as.matrix(cluster.averages[[assay]]@scale.data) %>% add.feature.names
  writeData(wb, sheet = "all_markers_scaled", all_markers_scaled, rowNames=TRUE)
  
  addWorksheet(wb, sheetName = "heatmap_markers_scaled")
  heatmap_markers_scaled <- as.matrix(cluster.averages[[assay]]@scale.data[markers.use$gene,]) %>% add.feature.names
  writeData(wb, sheet = "heatmap_markers_scaled", heatmap_markers_scaled, rowNames = TRUE)
  
  saveWorkbook(wb, file = xlsx.name, overwrite=TRUE)
  
  return(list(all_markers_lognorm = all_markers_lognorm, 
              heatmap_markers_lognorm = heatmap_markers_lognorm, 
              all_markers_scaled = all_markers_scaled, 
              heatmap_markers_scaled = heatmap_markers_scaled
  ))
  
}


create.and.export.clus.avg <- function(wnn.clusterings, ident_clustering_resolution, marker.set, 
                                       clustering_resolution, 
                                       assay){
 #assay = SCT if RNA
  cluster.averages <- create.cluster.avg(wnn.clusterings, ident_clustering_resolution)
  
  require(dplyr)
  min.lfg <- 0.0
  max.p.adj <- 0.05
  limit <- 20
  mask.mt.rp <- T
  markers.use=subset(marker.set,avg_log2FC > min.lfg & 
                       p_val_adj < max.p.adj & 
                       !(mask.mt.rp & (grepl("^RP[SL]", gene) | grepl("^MT", gene)))) %>%
    group_by(cluster) %>% 
    top_n(limit, -p_val_adj) %>% 
    top_n(limit, abs(pct.1-pct.2)) %>%
    arrange(as.factor(cluster))
  
  xlsx.name <- paste0("cluster.averages_", clustering_resolution, ".xlsx")
  list_clus_avg <- export.cluster.avg.and.heatmap(cluster.averages, markers.use, xlsx.name, assay = assay)
  #list_clus_avg$markers <- marker.set
  
  # fileName <- paste0("clus.avg", '_', clustering_resolution,'_', assay, ".pdf")
  # plot.heatmap.seurat(wnn.clusterings, marker.set, fileName)
  # 
  #list_clus_avg <- lapply(list_clus_avg, as.data.frame)
  
  return(list_clus_avg)
}

################

find.markers_wnn_0.6 <- function(seurat.wnn.clusterings){
  wnn_rna_0.6 <- find.markers.by.clustering(seurat.wnn.clusterings, "clusters_wsnn_0.6", assay="RNA")
  wnn_peaks_0.6 <- find.markers.by.clustering(seurat.wnn.clusterings, "clusters_wsnn_0.6", assay="peaks")
  return(list(rna = wnn_rna_0.6, 
              atac = wnn_peaks_0.6))
}


plot.heatmap_wnn_0.6_peaks_rna <- function(seurat.wnn.clusterings, wnn_0.6_markers){
  
  DefaultAssay(seurat.wnn.clusterings) <- "RNA"
  Idents(seurat.wnn.clusterings) <- seurat.wnn.clusterings@meta.data$clusters_wsnn_0.6
  fileName <- "heatmap.wnn_0.6.RNA.markers.pdf"
  clus.avg_wsnn_0.6_rna <- plot.heatmap.seurat(seurat.wnn.clusterings, 
                                               markers = wnn_0.6_markers$rna, 
                                               "RNA", fileName=fileName, marker.set.name="wnn_0.6_RNA")
  
  DefaultAssay(seurat.wnn.clusterings) <- "peaks"
  Idents(seurat.wnn.clusterings) <- seurat.wnn.clusterings@meta.data$clusters_wsnn_0.6
  fileName <- "heatmap.wnn_0.6.peaks.markers.pdf"
  clus.avg_wsnn_0.6_peaks <- plot.heatmap.seurat(seurat.wnn.clusterings, 
                                               markers = wnn_0.6_markers$atac, 
                                               "peaks", fileName=fileName, marker.set.name="wnn_0.6_peaks")
  
  
  return(list(rna = clus.avg_wsnn_0.6_rna, 
              atac = clus.avg_wsnn_0.6_peaks))
}


calculate.wnn_0.6 <- function(seurat.wnn.clusterings){
  
  wnn_0.6_markers <- find.markers_wnn_0.6(seurat.wnn.clusterings) 
  clus.avg_wsnn_0.6 <- plot.heatmap_wnn_0.6_peaks_rna(seurat.wnn.clusterings, wnn_0.6_markers)
  
  return(markers = wnn_0.6_markers, 
         clus.avg = clus.avg_wsnn_0.6)
}

# add.feature.names <- function(list_of_list){
#   update_lst <- lapply(list_of_list, function(dataset){
#     dataset <- dataset %>% as_tibble %>% add_column(feature.names=rownames(dataset), .before=0) %>%
#       as.data.frame()
#   })
#   return(update_lst)
# }

clus.avg.prep.for.h5 <- function(clus.avg_main.cat, marker.set){
  
  #add feature names
  modality <- c("rna", "atac")
  library(tibble)
  
  clus.avg_main.cat_h5 <- lapply(names(clus.avg_main.cat), function(modality){
    lst_resolutions <- clus.avg_main.cat[[modality]]
    lst_resolutions_update <- add.feature.names(lst_resolutions)
    return(lst_resolutions_update)
  })
  names(clus.avg_main.cat_h5) <- names(clus.avg_main.cat)
  
  #add marker list to this clus.avg
  for (layer in names(clus.avg_main.cat_h5)){
    clus.avg_main.cat_h5[[layer]][["markers"]] <- marker.set[[layer]]
  }
  return(clus.avg_main.cat_h5)
}



