# Update metadata for one field at a time in covid19_datasets.h5seurat
# usage: Rscript update_metadata.R <fieldname>
# eg: Rscript update_metadata.R disease_severity

library(SeuratDisk)
library(data.table)
setDTthreads(4)

# args = name of field to update
args = commandArgs(trailingOnly=TRUE)
field <- args[1]

h5seuratPath <- file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed", "covid19_datasets.h5seurat")
metadataPath <- file.path(Sys.getenv("DATA_DIR"), "metadata", "all_meta.tsv")

all_meta <- fread(metadataPath)
if (!field %in% names(all_meta)) stop("provide valid field to update")
hfile <- Connect(h5seuratPath, mode = "r+")
dataset_hfile <- hfile[["meta.data/dataset"]][]
sample_hfile <- hfile[["meta.data/sample"]][]
keys <- data.table(dataset = dataset_hfile, sample = sample_hfile)
new_values <- merge(keys, all_meta,  by = c("dataset", "sample"), all.x = TRUE, sort = FALSE)[[field]]
if (field %in% names(hfile[["meta.data"]])) {
  hfile[["meta.data"]][[field]][] <- new_values
} else {
  hfile[["meta.data"]][[field]] <- new_values
}
hfile$close_all()
