library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(Azimuth)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rhdf5)
library(dplyr)
library(SingleR)
library(future)

plan(sequential)
set.seed(1234)

resolution.lst <- c(0.05, 0.08, 0.6, 0.8, 1)
#c(0.05, 0.08, 0.6,0.8, 1, 1.2, 2, 4)

create_seurat_multiome <- function(data.path, meta.data=FALSE){
  #this Seurat is without meta.data
  #Updated 12/1/2025 to handle both H5 and Matrix Market (MTX) formats
  
  library(EnsDb.Hsapiens.v86)
  # get gene annotations for hg38
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  
  # Check for H5 file first
  h5_files <- list.files(data.path, pattern = ".*filtered_feature_bc_matrix\\.h5$", full.names = FALSE)
  
  if(length(h5_files) > 0){
    # === H5 Format ===
    path.10x <- paste0(data.path, "/", h5_files[1])
    print(paste("Using H5 file:", path.10x))
    counts <- Read10X_h5(path.10x)
    
    # Extract Gene Expression and Peaks
    rna_counts <- counts$`Gene Expression`
    atac_counts <- counts$Peaks
    
  } else {
    # === Matrix Market (MTX) Format ===
    # Look for MTX files
    mtx_files <- list.files(data.path, pattern = ".*\\.mtx(\\.gz)?$", full.names = FALSE)
    
    if(length(mtx_files) > 0){
      print("Using Matrix Market format (MTX + TSV files)")
      
      # Find barcodes and features files (handles custom prefixes like GSM5765322_100809M-TRM_)
      barcodes_file <- list.files(data.path, pattern = ".*barcodes\\.tsv(\\.gz)?$", full.names = TRUE)
      features_file <- list.files(data.path, pattern = ".*features\\.tsv(\\.gz)?$", full.names = TRUE)
      matrix_file <- list.files(data.path, pattern = ".*\\.mtx(\\.gz)?$", full.names = TRUE)
      
      if(length(barcodes_file) == 0 | length(features_file) == 0 | length(matrix_file) == 0){
        stop("MTX format requires: barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files")
      }
      
      # Manually read MTX files to handle custom prefixes
      library(Matrix)
      
      # Read barcodes
      barcodes <- readLines(gzfile(barcodes_file[1]))
      
      # Read features
      features <- read.table(gzfile(features_file[1]), sep="\t", stringsAsFactors=FALSE)
      
      # Read matrix (MTX format)
      mtx <- readMM(gzfile(matrix_file[1]))
      
      # Separate RNA and Peaks based on V3 column
      rna_idx <- which(features$V3 == "Gene Expression")
      peaks_idx <- which(features$V3 == "Peaks")
      
      if(length(rna_idx) > 0){
        rna_matrix <- mtx[rna_idx, ]
        rna_names <- features$V2[rna_idx]
        # Handle duplicate feature names
        if(sum(duplicated(rna_names)) > 0){
          print(paste("Warning: Found", sum(duplicated(rna_names)), "duplicate RNA feature names. Making unique."))
          rna_names <- make.unique(rna_names)
        }
        rownames(rna_matrix) <- rna_names
        colnames(rna_matrix) <- barcodes
      } else {
        stop("No Gene Expression features found in MTX data")
      }
      
      if(length(peaks_idx) > 0){
        peaks_matrix <- mtx[peaks_idx, ]
        peak_names <- features$V1[peaks_idx]
        # Handle duplicate feature names
        if(sum(duplicated(peak_names)) > 0){
          print(paste("Warning: Found", sum(duplicated(peak_names)), "duplicate Peak feature names. Making unique."))
          peak_names <- make.unique(peak_names)
        }
        rownames(peaks_matrix) <- peak_names
        colnames(peaks_matrix) <- barcodes
        atac_counts <- peaks_matrix
        print("Warning: No Peaks found in MTX data, ATAC assay will not be created")
        atac_counts <- NULL
      }
      
      rna_counts <- rna_matrix
      
    } else {
      # Fallback to directory format (filtered_feature_bc_matrix)
      path.10x <- paste0(data.path, "/filtered_feature_bc_matrix")
      print(paste("Using directory format:", path.10x))
      counts <- Read10X(path.10x)
      
      if(is.list(counts)){
        rna_counts <- counts$`Gene Expression`
        atac_counts <- counts$Peaks
      } else {
        rna_counts <- counts
        atac_counts <- NULL
      }
    }
  }
  
  # Create a Seurat object with RNA
  seurat <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA"
  )
  seurat@assays[["RNA"]]@layers[["counts"]]@Dimnames[[1]] <- rownames(rna_counts)
  seurat@assays[["RNA"]]@layers[["counts"]]@Dimnames[[2]] <- colnames(rna_counts)
  
  # Check for sample metadata file (created by merge_multiome_samples.R)
  metadata_file <- list.files(data.path, pattern = ".*sample_metadata\\.tsv$", full.names = TRUE)
  if(length(metadata_file) > 0){
    print(paste("Found sample metadata file:", basename(metadata_file[1])))
    sample_metadata <- read.delim(metadata_file[1], stringsAsFactors = FALSE)
    
    # Match barcodes to cells in Seurat object
    metadata_subset <- sample_metadata[match(colnames(seurat), sample_metadata$barcode), ]
    
    # Add to Seurat metadata
    if(nrow(metadata_subset) == ncol(seurat) && !all(is.na(metadata_subset$sample))){
      seurat$sample <- metadata_subset$sample
      seurat$tissue <- metadata_subset$tissue
      seurat$patient <- metadata_subset$patient
      print(paste("✓ Added sample metadata:", 
                  paste(unique(metadata_subset$sample[!is.na(metadata_subset$sample)]), collapse=", ")))
    } else {
      print("Warning: Sample metadata barcodes do not match Seurat object barcodes")
    }
  }
  
  # Add ATAC assay if available
  if(!is.null(atac_counts)){
    fragpath <- list.files(data.path, pattern='.*_fragments.tsv.gz', recursive = F)
    
    if(length(fragpath) > 0){
      # Create ATAC assay with fragments file
      seurat[["ATAC"]] <- CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = 'hg38',
        fragments = paste0(data.path,'/', fragpath[1]),
        annotation = annotation
      )
    } else {
      # Create ATAC assay without fragments file (still valid, just limited functionality)
      print("Warning: No fragments file found. ATAC assay created without fragment support.")
      seurat[["ATAC"]] <- CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = 'hg38',
        annotation = annotation
      )
    }
  } else {
    print("Note: No ATAC data found in this dataset")
  }
  
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
  seurat.raw <- CallPeaks(seurat.raw, 
                          outdir = data.path,
                          macs2.path = "/Users/chloetran/miniconda3/envs/macs3/bin/macs3")
  return(seurat.raw)
}

# Smart peak handling - tries CallPeaks, falls back to pre-called peaks if no fragments
call.peaks_smart <- function(seurat.raw, data.path){
  DefaultAssay(seurat.raw) <- "ATAC"
  
  # Try to call peaks with CallPeaks (requires fragments)
  tryCatch({
    cat("Attempting de novo peak calling with MACS3...\n")
    seurat.raw <- CallPeaks(seurat.raw, 
                            outdir = data.path,
                            macs2.path = "/Users/chloetran/miniconda3/envs/macs3/bin/macs3")
    cat("✓ De novo peak calling successful\n")
    return(seurat.raw)
  }, error = function(e) {
    cat("⚠ De novo peak calling failed:", conditionMessage(e), "\n")
    cat("Attempting to use pre-called peaks instead...\n")
    
    # Fallback: look for pre-called peak files
    peak_files <- list.files(data.path, pattern = ".*_peaks\\.bed(\\.gz)?$", full.names = TRUE)
    
    if(length(peak_files) > 0) {
      cat("Found pre-called peaks:", basename(peak_files[1]), "\n")
      peaks <- read.table(peak_files[1], sep="\t", stringsAsFactors=FALSE)
      colnames(peaks) <- c("chrom", "start", "end")
      cat("✓ Using", nrow(peaks), "pre-called peaks\n")
      
      # Return seurat object - caller should use create_peaks_assay()
      return(seurat.raw)
    } else {
      stop("No fragment files for peak calling AND no pre-called peak files found in:", data.path)
    }
  })
}

create_peaks_assay <- function(peaks, seurat){
  library(GenomicRanges)
  DefaultAssay(seurat) <- "ATAC"
  #input 'peaks' can be either from 'peak calling' from MACS2 or from 'peaks.bed' file
  peaks <- makeGRangesFromDataFrame(df = peaks)
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions

  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  # require(future)
  # plan("multisession", workers = 5)
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

# Use pre-calculated peak counts from MTX (no fragments needed)
use_precalculated_peaks <- function(seurat, peak_bed_path = NULL){
  library(GenomicRanges)
  DefaultAssay(seurat) <- "ATAC"
  
  # The ATAC assay already has peak counts from the MTX format
  # Peak names are in format "chr1-9963-10510" (chr-start-end)
  # We just need to validate/filter them
  
  cat("Using pre-calculated peak counts from ATAC assay\n")
  
  # Get current peak GRanges from row names
  peak_names <- rownames(seurat@assays$ATAC@counts)
  cat("Total peaks in ATAC assay:", length(peak_names), "\n")
  
  # Parse peak names (format: chr1-9963-10510)
  peak_df <- data.frame(
    peak = peak_names,
    chrom = sapply(strsplit(peak_names, "-"), "[", 1),
    start = as.numeric(sapply(strsplit(peak_names, "-"), "[", 2)),
    end = as.numeric(sapply(strsplit(peak_names, "-"), "[", 3)),
    stringsAsFactors = FALSE
  )
  
  # Convert to GRanges
  peaks_gr <- makeGRangesFromDataFrame(df = peak_df[, c("chrom", "start", "end")], keep.extra.columns = TRUE)
  names(peaks_gr) <- peak_df$peak  # Preserve original peak names
  
  # Remove peaks on nonstandard chromosomes
  peaks_gr <- keepStandardChromosomes(peaks_gr, pruning.mode = "coarse")
  
  # Remove blacklist regions
  peaks_gr <- subsetByOverlaps(x = peaks_gr, ranges = blacklist_hg38_unified, invert = TRUE)
  
  cat("After filtering (nonstandard chromosomes + blacklist):", length(peaks_gr), "peaks\n")
  
  # Subset ATAC counts to keep only these peaks (using names, not numeric indices)
  peaks_to_keep <- names(peaks_gr)
  atac_counts <- seurat@assays$ATAC@counts[peaks_to_keep, ]
  
  # Create new "peaks" assay with filtered peaks (mimics create_peaks_assay behavior)
  # Note: rownames should already be in correct format from atac_counts
  seurat[["peaks"]] <- CreateChromatinAssay(
    counts = atac_counts,
    annotation = Annotation(seurat)
  )
  
  cat("✓ Created 'peaks' assay with", length(peaks_gr), "high-quality peaks\n")
  return(seurat)
}

# Use ATAC assay directly as peaks (no filtering, no fragments needed)
# For merged samples where peaks are already union of both samples
use_atac_as_peaks <- function(seurat){
  library(GenomicRanges)
  DefaultAssay(seurat) <- "ATAC"
  
  cat("Using ATAC assay directly as 'peaks' (no filtering)\n")
  cat("Total peaks in ATAC assay:", nrow(seurat@assays$ATAC@counts), "\n")
  
  # Simply copy ATAC assay to peaks assay
  # This preserves all peaks including sample-specific ones from merged data
  seurat[["peaks"]] <- seurat[["ATAC"]]
  
  cat("✓ Created 'peaks' assay with", nrow(seurat@assays$peaks@counts), "peaks (identical to ATAC)\n")
  cat("Note: No filtering applied - all peaks retained for merged sample analysis\n")
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
    
    seurat <- call.peaks_from_signac(seurat.raw.ls[[data.name]], path)
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
  # require(SingleR)
  # require(future)
  DefaultAssay(seurat.raw) <- "RNA"
  raw.reads <- GetAssayData(seurat.raw, assay = "RNA", layer = "counts")
  #options(future.globals.maxSize = 800000 * 1024^2)
  preds <- sapply(names(refs), function(ref) {
    sapply(c("main","fine"), function(r) {
      SingleR(raw.reads, refs[[ref]], refs[[ref]][[paste0("label.",r)]] ) }
    )
  })
  return(preds)
}

#================
qc_ATAC_premerge <- function(obj) {
  DefaultAssay(obj) <- "ATAC"
  
  # Fragment-based QC metrics
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)
  
  # Reads/fragments in peaks
  obj$atac_peak_region_fragments <- Matrix::colSums(
    GetAssayData(obj, slot = "counts")
  )
  
  # Total fragments per cell
  frag_path <- Fragments(obj)[[1]]@path
  frag_df <- CountFragments(frag_path)
  rownames(frag_df) <- frag_df$CB
  
  obj$atac_fragments <- frag_df[colnames(obj), "frequency_count"]
  
  # Percent reads in peaks
  obj$pct_reads_in_peaks <-
    obj$atac_peak_region_fragments / obj$atac_fragments * 100
  
  # Filter
  obj <- subset(
    obj,
    subset =
      atac_peak_region_fragments > 3000 &
      atac_peak_region_fragments < 20000 &
      pct_reads_in_peaks > 15 &
      nucleosome_signal < 4 &
      TSS.enrichment > 2
  )
  
  return(obj)
}
#====================

# #############
#QC in ATAC
#QC in ATAC
qc_ATAC <- function(seurat, meta.data =TRUE){
  #meta.data means "atac_peak_region_fragments"
  #options(future.globals.maxSize = 30 * 1024^3)
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
  }
  if(meta.data==FALSE) { #if meta.data not avai, only QC ATAC based on nCount_ATAC
    seurat <- subset(
      x = seurat,
      subset = nCount_ATAC < 100000 &
        nCount_RNA < 25000 &
        nCount_ATAC > 1000 &
        nCount_RNA > 1000 &
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



# #################
integrate.seurat.by.RNA <- function(seurat.merge, nGenes, nCCA, split.by="patient"){
    # nCCA=30
  # nGenes=6000
  #library(BiocParallel)

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

  #future::plan("sequential")

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
  return(seurat)
}


do.clusterings <- function(seurat, assay, resolution.lst) {
  DefaultAssay(seurat) <- assay
  
  # Find available graphs
  available_graphs <- names(seurat@graphs)
  
  # Determine which graph to use - check what actually exists
  possible_names <- c(
    paste0(assay, "_snn"),
    paste0(assay, "_nn"),
    "SCT_snn",  # SCTransform creates SCT_snn
    "RNA_snn",
    "peaks_snn",
    "ATAC_snn"
  )
  
  graph.name <- NULL
  for(name in possible_names) {
    if(name %in% available_graphs) {
      graph.name <- name
      break
    }
  }
  
  if(is.null(graph.name)) {
    stop("No suitable graph found in Seurat object. Available graphs: ", 
         paste(available_graphs, collapse = ", "))
  }
  
  cat("Using graph:", graph.name, "\n")
  
  for(resolution in resolution.lst){
    seurat <- FindClusters(seurat, graph.name = graph.name, resolution = resolution)
    seurat[[paste0("clusters_",assay, "_",resolution)]] = seurat$seurat_clusters
  }
  return(seurat)
}




#############
#analyze ATAC


analyze_atac <- function(RWPE.combined, debatch=FALSE){

  DefaultAssay(RWPE.combined) <- "peaks"

  RWPE.combined <- RunTFIDF(RWPE.combined)
  RWPE.combined <- FindTopFeatures(RWPE.combined, min.cutoff = 20) #Only peaks presented in more than 20 cells were selected
  #to bypass the version conflicted of "Matrix" and "irlba"
  #RWPE.combined <- RunSVD(RWPE.combined)
  #this "svd.method" only avai with Matrix version 1.7-3
  RWPE.combined <- RunSVD(RWPE.combined, svd.method = "rsvd")
  
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
  # require(future)
  # options(future.globals.maxSize = 800000 * 1024^2)
  # plan("multisession", workers = 4)
  
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
    # 1) use all wsnn resolutions present in meta.data
    clusterings <- grep("clusters_wsnn_.*", colnames(seurat@meta.data), value = TRUE)
    
    # 2) append two extra "clusterings" (one row each in the final list)
    extra <- c(
      "curated.celltype",
      "reported.annotation",
      colnames(seurat@meta.data)[startsWith(colnames(seurat@meta.data), "azimuth")]
    )
    present_extra <- intersect(extra, colnames(seurat@meta.data))
    missing_extra <- setdiff(extra, present_extra)
    
    if (length(missing_extra) > 0) {
      warning(
        "Missing metadata columns (skipping): ",
        paste(missing_extra, collapse = ", ")
      )
    }
    
    clusterings <- c(clusterings, present_extra)
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
do.curated.cell.types <- function(cell.and.cluster.stats, default.clustering = "clusters_wsnn_0.8") {
  
  stopifnot(default.clustering %in% names(cell.and.cluster.stats))
  clustering <- default.clustering
  #clustering <- names(cell.and.cluster.stats)[1]
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
do.curated.cell.types2 <- function(
    cell.and.cluster.stats,
    mode = c("majority", "direct"),
    clustering = NULL,
    direct_label_col = NULL,
    cluster_col_for_suffix = NULL,
    exception_clusters = c(),
    append_cluster = TRUE
) {
  mode <- match.arg(mode)
  
  # -------------------------
  # MAJORITY MODE
  # -------------------------
  if (mode == "majority") {
    
    stopifnot(!is.null(clustering))
    stopifnot(clustering %in% names(cell.and.cluster.stats))
    
    cell_stats <- cell.and.cluster.stats[[clustering]]$cell.stats
    cluster_stats <- cell.and.cluster.stats[[clustering]]$cluster.stats
    
    major_ct <- as.character(cluster_stats$ct.main)
    names(major_ct) <- as.character(cluster_stats$clusterN)
    
    cluster_vec <- cell_stats[[clustering]]
    cluster_chr <- as.character(cluster_vec)
    
    ct <- as.character(major_ct[cluster_chr])
    ct[cluster_vec %in% exception_clusters] <- NA_character_
    
    if (append_cluster) {
      return(paste0(ct, "-", cluster_chr))
    } else {
      return(ct)
    }
  }
  
  # -------------------------
  # DIRECT MODE
  # -------------------------
  if (mode == "direct") {
    
    stopifnot(!is.null(direct_label_col))
    
    # find the direct label inside any clustering
    found <- FALSE
    for (cl in names(cell.and.cluster.stats)) {
      cs <- cell.and.cluster.stats[[cl]]$cell.stats
      if (!is.null(colnames(cs)) && direct_label_col %in% colnames(cs)) {
        direct_ct <- as.character(cs[[direct_label_col]])
        found <- TRUE
        break
      }
    }
    
    if (!found) stop("direct_label_col not found in any cell.stats")
    
    direct_ct[direct_ct == ""] <- NA_character_
    
    if (append_cluster) {
      stopifnot(!is.null(cluster_col_for_suffix))
      
      # find cluster column for suffix
      found_cluster <- FALSE
      for (cl in names(cell.and.cluster.stats)) {
        cs <- cell.and.cluster.stats[[cl]]$cell.stats
        if (cluster_col_for_suffix %in% colnames(cs)) {
          cluster_vec <- as.character(cs[[cluster_col_for_suffix]])
          found_cluster <- TRUE
          break
        }
      }
      
      if (!found_cluster) stop("cluster_col_for_suffix not found")
      
      return(paste0(direct_ct, "-", cluster_vec))
    }
    
    return(direct_ct)
  }
}
##########################################

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
  
  # Add consistent metadata from seurat@meta.data
  # These columns should match what's in match_curated_celltyes_multiome_discovery.R
  
  # sample.type (from orig.ident if not already present)
  if ("sample.type" %in% colnames(seurat@meta.data)) {
    covs$sample.type <- seurat@meta.data$sample.type
  } else if ("orig.ident" %in% colnames(seurat@meta.data)) {
    covs$sample.type <- seurat@meta.data$orig.ident
  }
  
  # patient ID
  if ("patient" %in% colnames(seurat@meta.data)) {
    covs$patient <- seurat@meta.data$patient
  }
  
  # sample (full sample name like "100809M-ReCir")
  if ("sample" %in% colnames(seurat@meta.data)) {
    covs$sample <- seurat@meta.data$sample
  }
  
  # curated.celltype (curated cell type annotations)
  if ("curated.celltype" %in% colnames(seurat@meta.data)) {
    covs$curated.celltype <- seurat@meta.data$curated.celltype
  }
  
  # reported.annotation
  if ("reported.annotation" %in% colnames(seurat@meta.data)) {
    covs$reported.annotation <- seurat@meta.data$reported.annotation
  }
  # azimuth_predicted.celltype.l2
  az_cols <- grep("^azimuth", colnames(seurat@meta.data), value = TRUE)
  
  if (length(az_cols) > 0) {
    covs[az_cols] <- seurat@meta.data[, az_cols, drop = FALSE]
  }
  return(covs)
}


################

get.cell.and.cluster.stats.for.clustering <- function(seurat, preds, clustering="clusters_wsnn_0.8") {
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
  meta_cols <- colnames(seurat@meta.data)
  clusterings <- meta_cols[grep("^clusters_|^azimuth", meta_cols)]

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
annotate.curated.cell.types_2assays <- function(curated.cell.types, markers_sct_clusterings_update) {
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
  
  default.markers.lst <- markers_sct_clusterings_update[["default.cts"]]
  
  order_default.markers.lst <- lapply(default.markers.lst, function(markers){
    markers$cluster <- factor(markers$cluster, levels = cts.order)
    markers <- markers %>% arrange(cluster)
    
  })
  
  return(list(cts.order=cts.order, 
              markers=order_default.markers.lst))
}



###################
default.cts_find.markers.and.add.to.markers.sct.clusterings <- function(seurat.wnn.clusterings, default.cts, markers_sct_clusterings){
  seurat.wnn.clusterings$default.cts <- default.cts
  markers_default.cts_2assays <- find.markers.for.two.assays(seurat.wnn.clusterings, "default.cts", resolution.lst)
  save(markers_default.cts_2assays, file = "markers_default.cts_2assays.RData")
  
  markers_sct_clusterings$default.cts <- markers_default.cts_2assays
  return(markers_sct_clusterings)
}
##############

###############

default.cts_find.markers.and.add <- function(
    seurat.wnn.clusterings,
    default.cts,
    markers_lists,           # named list, e.g. list(sct=markers_sct_clusterings, peaks=markers_peaks_clusterings)
    resolution.lst,
    save_file = NULL         # e.g. "markers_default.cts_2assays.RData" or NULL to skip saving
) {
  # Add default.cts to Seurat meta
  seurat.wnn.clusterings$default.cts <- default.cts
  
  # Compute markers once (should return list(rna=..., atac=...) or similar)
  markers_default.cts_2assays <- find.markers.for.two.assays(
    seurat.wnn.clusterings,
    "default.cts",
    resolution.lst
  )
  
  # Optional save
  if (!is.null(save_file)) {
    save(markers_default.cts_2assays, file = save_file)
  }
  
  # Update each marker container
  for (nm in names(markers_lists)) {
    markers_lists[[nm]]$default.cts <- markers_default.cts_2assays
  }
  
  markers_lists
}

##############
###############

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
write_mm_h5 <- function(seurat, covs_discrete, heatmap_artifacts, filename,join_rna) {
  
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
  write_h5(filename, seurat, layers = layers, assays = c("RNA", "peaks", "GeneActivity"),join_rna)
           
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
    # Extract base cell type name (before the hyphen and cluster number)
    # e.g., "B cells-1" -> "B cells", "Monocytes-2" -> "Monocytes"
    base_ct <- sub("-[0-9]+$", "", cts)
    
    sub.categories <- dplyr::case_when(
      grepl("^B[- ].*", base_ct) ~ "B-cells", 
      
      grepl("^CD8\\+", base_ct) ~ "CD8",
      grepl("^CD4\\+", base_ct) ~ "CD4",
      
      grepl("CD14", base_ct) ~ "Monocytes",
      grepl("CD16", base_ct) ~ "Monocytes",
      grepl("Monocytes", base_ct) ~ "Monocytes",
      grepl("ASDC", base_ct) ~ "Myeloids", 
      grepl("Macrophages", base_ct) ~ "Myeloids",
      
      grepl("HSPC", base_ct) ~ "Others",
      grepl("Eryth", base_ct) ~ "Others",
      grepl("ILC", base_ct) ~ "Others",
      
      grepl("^NK", base_ct) ~ "NK cells",
      
      TRUE ~ cts  # Keep original if no match
    )
    return(sub.categories)
    }

do.main.categories <- function(sub.categories){
  main.categories <- dplyr::case_when(
    grepl("^B[- ]", sub.categories) ~ "Lymphoids", 
    grepl("^NK", sub.categories) ~ "Lymphoids",
    
    grepl("gdT", sub.categories) ~ "Lymphoids",
    grepl("dnT", sub.categories) ~ "Lymphoids",
    
    grepl("CD8", sub.categories) ~ "Lymphoids",
    grepl("CD4", sub.categories) ~ "Lymphoids",
    grepl("Treg", sub.categories) ~ "Lymphoids",
    
    grepl("Monocytes", sub.categories) ~ "Monocytes", 
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
  
  return(list(seurat.wnn = seurat.wnn,
              main.categories = list(main.categories = main.categories, 
                                     main.cat_rna.markers = rna.main.cat, 
                                     main.cat_peak.markers = peaks.main.cat
                                     ), 
              sub.categories = list(sub.categories = sub.categories, 
                                    sub.cat_rna.markers = rna.sub.cat, 
                                    sub.cat_peak.markers = peaks.sub.cat)))
}



# azimuth annotation
# azimuth annotation
azimuth_annotation <- function(
    obj,
    reference,
    pred_cols = NULL,
    assay = "RNA",
    join_layers = TRUE,
    verbose = TRUE
) {
  stopifnot(inherits(obj, "Seurat"))
  
  # reference validation (handle missing() safely)
  if (missing(reference)) stop("`reference` is required.")
  if (is.null(reference)) stop("`reference` is required and cannot be NULL.")
  if (is.character(reference) && !nzchar(reference)) stop("`reference` cannot be empty.")
  
  # assay check
  if (!assay %in% names(obj@assays)) {
    stop(sprintf(
      "Assay '%s' not found. Available assays: %s",
      assay, paste(names(obj@assays), collapse = ", ")
    ))
  }
  
  # RunAzimuth availability
  if (!exists("RunAzimuth", mode = "function")) {
    stop("RunAzimuth() not found. Try: library(Azimuth)")
  }
  
  DefaultAssay(obj) <- assay
  
  # Join layers (Seurat v5)
  if (isTRUE(join_layers)) {
    obj[[assay]] <- tryCatch({
      if (exists("Layers", mode = "function")) {
        lyr <- Layers(obj[[assay]])
        if (length(lyr) > 1) {
          if (isTRUE(verbose)) {
            message("Joining layers in assay '", assay, "': ", paste(lyr, collapse = ", "))
          }
          JoinLayers(obj[[assay]])
        } else {
          obj[[assay]]
        }
      } else {
        # If Layers() doesn't exist, still attempt JoinLayers
        JoinLayers(obj[[assay]])
      }
    }, error = function(e) {
      warning("JoinLayers() failed; continuing without joining layers. Error: ", conditionMessage(e))
      obj[[assay]]
    })
  }
  
  # a safe tag used only for naming output columns
  ref_tag <- if (is.character(reference)) reference else class(reference)[1]
  ref_tag <- gsub("[^A-Za-z0-9_]+", "_", ref_tag)
  
  if (isTRUE(verbose)) message("Running Azimuth mapping with reference = ", ref_tag, " ...")
  obj <- RunAzimuth(query = obj, reference = reference)
  
  md <- colnames(obj[[]])
  
  # auto-select pred_cols if not provided
  if (is.null(pred_cols)) {
    if (all(c("predicted.class", "predicted.subclass", "predicted.cluster") %in% md)) {
      pred_cols <- c(
        "predicted.class",
        "predicted.subclass",
        "predicted.cluster",
        "predicted.cross_species_cluster"
      )
      pred_cols <- pred_cols[pred_cols %in% md]
    } else if (any(c("predicted.celltype.l1", "predicted.celltype.l2", "predicted.celltype.l3") %in% md)) {
      pred_cols <- c("predicted.celltype.l1", "predicted.celltype.l2", "predicted.celltype.l3")
      pred_cols <- pred_cols[pred_cols %in% md]
    } else if (any(c("predicted.annotation.l1", "predicted.annotation.l2", "predicted.annotation.l3") %in% md)) {
      pred_cols <- c("predicted.annotation.l1", "predicted.annotation.l2", "predicted.annotation.l3")
      pred_cols <- pred_cols[pred_cols %in% md]
    } else {
      # fallback: grab anything starting with "predicted." or "annotation."
      pred_cols <- grep("^(predicted\\.|annotation\\.)", md, value = TRUE)
    }
  }
  
  if (length(pred_cols) == 0) {
    stop("No prediction columns selected. Check `pred_cols` or inspect colnames(obj[[]]).")
  }
  
  missing_cols <- setdiff(pred_cols, md)
  if (length(missing_cols) > 0) {
    stop("These pred_cols are missing from metadata: ", paste(missing_cols, collapse = ", "))
  }
  
  # add columns
  for (pc in pred_cols) {
    out_col <- paste0("azimuth_", ref_tag, "_", pc)
    vals <- obj[[pc]][, 1]
    
    if (length(vals) != ncol(obj)) stop("Length mismatch for column: ", pc)
    
    if (out_col %in% md) {
      warning("Metadata column '", out_col, "' already exists and will be overwritten.")
    }
    
    obj[[out_col]] <- vals
    if (isTRUE(verbose)) message("Added metadata column: ", out_col)
  }
  
  obj
}




################################



