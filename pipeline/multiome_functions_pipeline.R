library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)

library(BSgenome.Hsapiens.UCSC.hg38)
library(rhdf5)
library(dplyr)
library(SingleR)
library(future)


set.seed(1234)
resolution.lst <- c(0.05, 0.08, 0.6,0.8)
#c(0.05, 0.08, 0.6,0.8, 1, 1.2, 2, 4)

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
  seurat@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]] <- counts[["Gene Expression"]]@Dimnames[[1]]
  seurat@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]] <- counts[["Gene Expression"]]@Dimnames[[2]]
  
  
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

call.peaks_from_signac<- function(seurat.raw, data.path){
  DefaultAssay(seurat.raw) <- "ATAC"
  seurat.raw <- CallPeaks(seurat.raw, outdir = data.path)
  return(seurat.raw)
}

create_peaks_assay <- function(peaks, seurat){
  library(GenomicRanges)
  DefaultAssay(seurat) <- "ATAC"
  #input 'peaks' can be either from 'peak calling' from MACS2 or from 'peaks.bed' file
  peaks <- makeGRangesFromDataFrame(df = peaks)
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions

  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  require(future)
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
    fragments = Fragments(seurat),
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
do.singler <- function(seurat.raw, refs) {
  require(SingleR)
  require(future)
  DefaultAssay(seurat.raw) <- "RNA"
  raw.reads <- GetAssayData(seurat.raw)
  options(future.globals.maxSize = 800000 * 1024^2)
  preds <- sapply(names(refs), function(ref) {
    sapply(c("main","fine"), function(r) {
      SingleR(raw.reads, refs[[ref]], refs[[ref]][[paste0("label.",r)]] ) }
    )
  })
  return(preds)
}



#############
#QC in ATAC
qc_ATAC <- function(seurat, meta.data =TRUE){
  #meta.data means "atac_peak_region_fragments"
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
  seurat <- seurat[,seurat[["percent.mt"]] <= 20 & seurat[["nFeature_RNA"]] >= 500]
  return(seurat)
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


do.clusterings <- function(seurat, assay, resolution.lst) {
  DefaultAssay(seurat) <- assay
 
 for(resolution in resolution.lst){
    seurat <- FindClusters(seurat, resolution = resolution)
    seurat[[paste0("clusters_",assay, "_",resolution)]] = seurat$seurat_clusters
  }
  return(seurat)
}


do.clusterings.WNN <- function(seurat, resolution.lst) {
  #no need to specify assay since WNN is already combining rna and atac
  #DefaultAssay(seurat) <- assay
  
    for(resolution in resolution.lst ) {
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
  #modality.weight.name = "RNA.weight",
  modality.weight.name = c("SCT.weight", "peaks.weight"),
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
  #https://github.com/satijalab/seurat/issues/1501
  #https://github.com/hbctraining/scRNA-seq_online/issues/14
  if ((assay=="SCT") || (assay =="integrated")){
    #when doing clustering, clustering based on SCT, which store < 36k genes 
    #but when finding marker genes, maybe should use RNA assay which store 36k genes?
    #but RNA assay doesnt have normalized counts like SCT
    
    #INPUT assay here is "SCT" because using SCT for clustering, thus clustering column name included "SCT" word
    assay <- "SCT"
    DefaultAssay(seurat) <- "SCT"
  }else{
    DefaultAssay(seurat) <- assay #in case of "peaks"
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

find.markers.for.two.assays <- function(seurat.wnn.clusterings, clustering, resolution.lst){
  rna.markers <- find.markers.by.clustering(seurat.wnn.clusterings, clustering, assay="SCT")
  peaks.markers <- find.markers.by.clustering(seurat.wnn.clusterings, clustering, assay="peaks")
  return(list(rna = rna.markers, 
              atac = peaks.markers))
}

find.markers.by.all.clusterings.for.two.assays <- function(seurat, group=0, assay, resolution.lst) {
  
  mid.point <- length(resolution.lst)/2
  if (group ==1){
    #first half of resolution list
    clusterings <- sapply(resolution.lst[1:mid.point], function(res){
              paste0("clusters_", assay, "_",res)
    })
  }
  if (group ==2){
    #second half of resolution list
    clusterings <- sapply(resolution.lst[(mid.point+1):length(resolution.lst)], function(res){
          paste0("clusters_", assay, "_",res)
})
  }
  if (group ==0){
    clusterings <- sapply(resolution.lst, function(res){
          paste0("clusters_", assay, "_",res)
      })
   
  }
  if(assay=="wnn"){
    clusterings <- grep("clusters_wsnn_.*", colnames(seurat@meta.data), value = TRUE)
  }

  ret <- lapply(clusterings, function(clustering) {
    print(paste0('start with ', clustering))
   # find.markers.by.clustering(seurat,clustering, assay)
    
    find.markers.for.two.assays(seurat, clustering)
    
  })
  names(ret) <- clusterings
  ret
}

formating_markers_clusterings <- function(markers_sct_clusterings) {
  ### in marker df, convert cluster column from factor to character
  #so when writing h5 file, these cluster names will be conserved (instead of being converted to integer by h5)
  markers_sct_clusterings_test <- lapply(names(markers_sct_clusterings), function(clustering){
    
    clustering_ret <- lapply(names(markers_sct_clusterings[[clustering]]), function(modality){
      print(clustering)
      print(modality)
      markers_df <- markers_sct_clusterings[[clustering]][[modality]]
      markers_df$cluster <- as.character(markers_df$cluster)
      markers_sct_clusterings[[clustering]][[modality]] <- markers_df
      return(markers_sct_clusterings[[clustering]][[modality]])
    })
    names(clustering_ret) <-  c("rna", "atac")
    return(clustering_ret)
    
  })
  names(markers_sct_clusterings_test) <- names(markers_sct_clusterings)
  return(markers_sct_clusterings_test)
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

#######################

###############################
do.curated.cell.types <- function(cell.and.cluster.stats) {
  
  clustering <- names(cell.and.cluster.stats)[1]
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


###################
default.cts_find.markers.and.add.to.markers.sct.clusterings <- function(seurat.wnn.clusterings, default.cts, markers_sct_clusterings){
  seurat.wnn.clusterings$default.cts <- default.cts
  markers_default.cts_2assays <- find.markers.for.two.assays(seurat.wnn.clusterings, "default.cts")
  save(markers_default.cts_2assays, file = "markers_default.cts_2assays.RData")
  
  markers_sct_clusterings$default.cts <- markers_default.cts_2assays
  return(markers_sct_clusterings)
}
##############

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

######################

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



create.heatmap.artifacts.for.clustering.in.two.modalities <- function(seurat.wnn.clusterings, markers_sct_clusterings, clustering){
  #for each clustering, create cluster_average in 2 modalities (rna and atac)
  artifacts_markers_sct_0.6_rna <- heatmap_artifacts_from_seurat_V5(seurat=seurat.wnn.clusterings, 
                                                                    clustering = clustering, 
                                                                    markers = markers_sct_clusterings[[clustering]]$rna, 
                                                                    assay="SCT")
  artifacts_markers_sct_0.6_peaks <- heatmap_artifacts_from_seurat_V5(seurat=seurat.wnn.clusterings, 
                                                                      clustering = clustering, 
                                                                      markers = markers_sct_clusterings[[clustering]]$atac, 
                                                                      assay="peaks")
  return(list(rna = artifacts_markers_sct_0.6_rna, 
              atac = artifacts_markers_sct_0.6_peaks))
}


create.artifacts.ALL.clusterings <- function(seurat, all.markers){
  ##create artifacts for ALL clusterings
  clusterings_artifacts <- lapply(names(all.markers), function(clustering){
    print(clustering)
    #check if Idents has more than 1 level
    if (length(names(table(seurat[[clustering]]))) >= 2){
      artifacts <- create.heatmap.artifacts.for.clustering.in.two.modalities(seurat = seurat, 
                                                                             markers = all.markers, 
                                                                             clustering = clustering)
    }
  })
  names(clusterings_artifacts) <- names(all.markers)
  #remove empty list
  clusterings_artifacts <- clusterings_artifacts[unlist(lapply(clusterings_artifacts, length) != 0)]
  return(clusterings_artifacts)
}

###########
assemble_heatmap_artifacts <- function(clus.avg.for.sct.clusterings, 
                                       clus.avg.for.peaks.clusterings, 
                                       clus.avg.for.wnn.clusterings){
  
  return(list(rna =clus.avg.for.sct.clusterings,
              atac = clus.avg.for.peaks.clusterings,
              wnn = clus.avg.for.wnn.clusterings))
}
############

##################
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
#######################
calculate_gene_activity <- function(seurat) {
  # # Set annotation for peaks in the ATAC assay
  # Annotation(seurat@assays$peaks) <- seurat@assays$ATAC@annotation
  # 
  # library(EnsDb.Hsapiens.v86)
  # # get gene annotations for hg38
  # annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  # seqlevelsStyle(annotation) <- "UCSC"
  
  
  # Compute gene activities
  gene_activities <- GeneActivity(seurat, assay = "peaks")
  
  # Add gene activities to Seurat object
  seurat[["GeneActivity"]] <- CreateAssayObject(counts = gene_activities)
  
  return(seurat)
}


attach_closest_features_to_atac_markers <- function(seurat, heatmap_artifacts) {
  # # Setting annotation for ATAC assay
  # Annotation(seurat@assays$peaks) <- seurat@assays$ATAC@annotation
  
  # Calculate closest feature
  closest_feature <- ClosestFeature(seurat@assays$peaks, regions = granges(seurat@assays$peaks))
  
  # Annotate markers function
  annotate_markers <- function(markers) {
    markers %>%
      left_join(closest_feature[, c("query_region", "gene_name", "type", "distance")],
                by = c("gene" = "query_region")) %>%
      select(gene, gene_name, type, distance, everything()) %>%
      mutate_if(is.factor, as.character)
  }
  
  # Iterate over each layer and clustering to annotate markers
  for (layer in names(heatmap_artifacts)) {
    print(paste0('working on ', layer))
    for (clustering in names(heatmap_artifacts[[layer]]$clusterings)) {
      current_clustering <- heatmap_artifacts[[layer]]$clusterings[[clustering]]
      print(paste0('   clustering: ', clustering))
      # If markers exist in the current clustering and layer is ATAC, annotate markers
      if ("markers" %in% names(current_clustering) && layer == "atac") {
        print('annotating single modal atac markers')
        heatmap_artifacts[[layer]]$clusterings[[clustering]]$markers <- annotate_markers(current_clustering$markers)
      } 
      
      # If ATAC is present in the current clustering, annotate markers
      else if ("atac" %in% names(current_clustering)) {
        heatmap_artifacts[[layer]]$clusterings[[clustering]]$atac$markers <- annotate_markers(current_clustering$atac$markers)
      }
    }
  }
  return(heatmap_artifacts)
}


#############
write_mm_h5 <- function(seurat, covs_discrete, heatmap_artifacts, filename) {
  
  # Calculate gene activity and attach closest features
  seurat <- seurat %>% calculate_gene_activity()
  heatmap_artifacts <- attach_closest_features_to_atac_markers(seurat, heatmap_artifacts)
  
  # Define QC features and clustering features
  qc_features <- c("nFeature_RNA", "nCount_RNA", "nFeature_ATAC", "nCount_ATAC",
                   "nFeature_peaks", "nCount_peaks", "percent.mt", "percent.rp",
                   "nucleosome_signal", "nucleosome_percentile", "TSS.enrichment", "TSS.percentile")
  
  clustering_features <- grep("clusters.*", colnames(seurat@meta.data), value = TRUE)
  
  # Define reductions
  reductions <- list(rna = c("tsne.rna", "umap.rna"),
                     atac = c("tsne.atac", "umap.atac"),
                     wnn = c("wnn.tsne", "wnn.umap"))
  
  # Process each layer and attach heatmap artifacts
  layer_names <- c("rna", "atac", "wnn")
  layers <- lapply(layer_names, function(layer) {
    layer_artifacts <- layer_artifacts_from_seurat(
      seurat, 
      qc_features = qc_features, 
      umap = reductions[[layer]][2],
      tsne = reductions[[layer]][1],
      clustering_names = clustering_features,
      extra_discrete_covs = covs_discrete
    )
    layer_artifacts[["clusterings"]] <- heatmap_artifacts[[layer]]
    layer_artifacts
  })
  
  names(layers) <- layer_names
  
  # Write to .h5 file
  write_h5(filename, seurat, layers = layers, assays = c("RNA", "peaks", "GeneActivity"))
           
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
      grepl("Macrophages", cts) ~ "Myeloids",
      
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



