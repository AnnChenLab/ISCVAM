start_time <- Sys.time()
print(paste0("starting time: ", start_time))


# library(future)
library(here)
library(tibble)

# Define base paths once

pipeline.dir <- here("pipelines", "sc_rna_pipeline")
code.dir <- file.path(pipeline.dir, "code")
ref.dir <- file.path(pipeline.dir, "reference_data")
results.dir <- file.path(pipeline.dir, "example.results")
data.dir <- file.path(pipeline.dir, "example.datasets")
# Use environment variable for data path, with fallback default
data.path <- data.dir

# Load sources and references
source(file.path(code.dir, "scRNA_functions_pipeline.R"))
source(file.path(code.dir, "iscvamR_h5_customized.R"))
load(file.path(ref.dir, "singler.refs.rdata"))

# options(future.globals.maxSize = 65000 * 1024^2)

meta.tisch <- read.csv(file.path(ref.dir, "meta.data_check.endo.track.progress.man.csv"))

#######
#array input from slurm, these will be dataset index
args <- commandArgs(trailingOnly = TRUE)
#input <- args[1]
input <- 1

dataset <- meta.tisch[[input, "Dataset.Name"]]


print(paste0("dataset number: ", input))
print(paste0("dataset name: ", dataset))
print(paste0("number of cells: ", meta.tisch[[input, "Cells"]]))

print(paste0("current directory: ", getwd()))


mainDir <- results.dir

#create a folder for each dataset (recursive = TRUE creates parent dirs if needed)
if (!dir.exists(file.path(mainDir, dataset))) {
  dir.create(file.path(mainDir, dataset), recursive = TRUE)
}


# Set the custom storage directory
setwd(file.path(mainDir, dataset))
print(paste0("project folder: ", getwd()))

fn.h5 <- paste0(dataset, ".test.pipeline.h5")
#
# #####
# #run 1st stage analysis for all cells
# Otherwise, create from TISCH inputs and save for reuse.

seurat.raw <- create.seurat.from.Tisch.h5(data.path, dataset)
save(seurat.raw, file = "seurat.raw.RData")

seurat.raw.clean <- remove_umap_from_tisch(seurat.raw)
save(seurat.raw.clean, file = "seurat.raw.clean.RData")

singler <- do.singler(seurat.raw.clean@assays[["RNA"]]@layers[["counts"]], refs)
save(singler, file = "singler.RData")

seurat.qc <- qc.seurat(seurat.raw.clean, "Human")
save(seurat.qc, file = "seurat.qc.RData")

seurat.sct <- seurat.SCT(seurat.qc, debatch=FALSE)
save(seurat.sct, file = "seurat.sct.RData")

seurat.pca <- seurat.PCA(seurat.sct, pcs.to.compute = 200)
save(seurat.pca, file = "seurat.pca.RData")

seurat.neighbors <- seurat.findNeighbors(seurat.pca, pcs.to.analyze=40)
save(seurat.neighbors, file = "seurat.neighbors.RData")

seurat.clusters <- seurat.findClusters(seurat.neighbors)
save(seurat.clusters, file = "seurat.clusters.RData")

seurat.viz <- seurat.visualization(seurat.clusters)
save(seurat.viz, file = "seurat.viz.RData")

seurat.clusterings <- do.clusterings(seurat.viz, "SCT")
save(seurat.clusterings, file = "seurat.clusterings.RData")

all.markers <- find.markers.by.all.clusterings(seurat.clusterings)
save(all.markers, file = "all.markers.RData")

export.rna.markers <- write_all.markers("all.markers.xlsx", all.markers)
save(export.rna.markers, file = "export.rna.markers.RData")

#load("/Users/chloetran/Documents/iscvam/pipeline/tisch.pipeline/res_local_running/acral/seurat.clusterings.RData")
tisch.markers <- find.markers.Tisch.major.cts(seurat.clusterings, dataset)
save(tisch.markers, file = "tisch.markers.RData")

cell.and.cluster.stats <- get.cell.and.cluster.stats.for.all.clusterings(seurat.clusterings, singler)
save(cell.and.cluster.stats, file = "cell.and.cluster.stats.RData")

export.cells.clusters.stats <- write_cells.clusters.stats('cells.and.clusterings.xlsx', cell.and.cluster.stats)
save(export.cells.clusters.stats, file = "export.cells.clusters.stats.RData")

curated.cell.types <- do.curated.cell.types(cell.and.cluster.stats)
save(curated.cell.types, file = "curated.cell.types.RData")

covs.clinical <- get.Tisch.covs.clinical(data.path, dataset, seurat.clusterings)
save(covs.clinical, file = "covs.clinical.RData")

covs <- get.covs(seurat.clusterings, singler, "Human", covs.clinical, curated.cell.types)
save(covs, file = "covs.RData")

seurat.list <- create.seurat.list.for.ISCVAM(seurat.clusterings, covs)

save(seurat.list, file = "seurat.list.RData")

##########################################################################################
#########################2nd stage analysis starting from here############################
##########################################################################################

seurat.list <- update.seurat.list(seurat.list, all.markers, tisch.markers)

seurat.list.all.layers <- do.multiStages(seurat.list)
save(seurat.list.all.layers, file = "seurat.list.multi.stages.RData")

layers <- create.artifacts.for.all.layers(seurat.list.all.layers)
save(layers, file = "all.layers.RData")

print(" ")
print("Start writing h5 file:")
write.iscvam.h5 <- writing.iscvam.h5(seurat.list.all.layers, layers, dataset, fn.h5)

end_time <- Sys.time()

running_time <- end_time - start_time
print(paste0("running time for this script: ", running_time))

