message(">>> merge_sce.R")
library(SummarizedExperiment)
library(Seurat)
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
library(SeuratDisk)

h5seuratDir <- file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed")
h5SCEDir <- file.path(h5seuratDir, "h5sce")

# 1. Load objects to merge into memeory as Seurat objects
message(">>> loading bigSeu...")
bigSeu <- LoadH5Seurat(file.path(h5seuratDir, "bigSeu_tomerge.h5seurat"))
message(">>> loading stephenson...")
stephenson <- LoadH5Seurat(file.path(h5seuratDir, "stephenson_2021_tomerge.h5seurat"))

# 2. Convert to SCE objects and save as HDF5
  # Need to add ADT separately as AltExp
message(">>> converting to SCE...")
bigSce <- as.SingleCellExperiment(bigSeu)
altExp(bigSce, "predicted_ADT") <- SummarizedExperiment(assays = list(logcounts = bigSeu[["predicted_ADT"]]@data))

steSce <- as.SingleCellExperiment(stephenson)
altExp(steSce, "predicted_ADT") <- SummarizedExperiment(assays = list(logcounts = stephenson[["predicted_ADT"]]@data))

message(">>> saving HDF5SummarizedExperiments...")
# extract_sparse_array on DelayedAbind with sparse backends is not yet implemented
# https://www.gitmemory.com/issue/Bioconductor/DelayedArray/80/743468429
bigSce <- saveHDF5SummarizedExperiment(bigSce, h5SCEDir, prefix = "big_", verbose = TRUE, as.sparse = FALSE, replace = FALSE)
steSce <- saveHDF5SummarizedExperiment(steSce, h5SCEDir, prefix = "ste_", verbose = TRUE, as.sparse = FALSE, replace = FALSE)

message(">>> combining...")
# Combine and save as hdf5 array
bigSce$nCount_predicted_ADT <- NULL
bigSce$nFeature_predicted_ADT <- NULL
combined <- cbind(bigSce, steSce)


# 4. Save the result
message(">>> Saving combined to ", h5SCEDir)
saveHDF5SummarizedExperiment(combined,
                             h5SCEDir,
                             as.sparse = TRUE,
                             prefix = "combined_",
                             replace = TRUE,
                             verbose = TRUE)

