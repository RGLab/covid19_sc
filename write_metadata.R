# Extract metadata from h5seurat file and write as tsv
# Ex: Rscript write_metadata.R liao_2020
library(SeuratDisk)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

dataset <- args[1]
datadir <- "/fh/fast/gottardo_r/ytian_working/covid19_datasets/h5seurat"

# Add sample field with timepoint to silvin dataset
if (dataset == "silvin_2020") {
  hfile <- Connect(file.path(datadir, paste0(dataset, "_processed.h5seurat")), mode = "r+")
  if (!"sample" %in% names(hfile[["meta.data"]])) {
    patient <- hfile[["meta.data/Characteristics.individual."]][]
    timepoint <- hfile[["meta.data/Factor.Value.sampling.time.point."]][]
    sample <- paste0(trimws(patient), "_", timepoint)
    hfile[["meta.data/sample"]] <- sample
  }
  hfile$close_all()
}


hfile <- Connect(file.path(datadir, paste0(dataset, "_processed.h5seurat")))
# https://hhoeflin.github.io/hdf5r/
metadataInfo <- hfile[["meta.data"]]$ls()
metadataList <- lapply(metadataInfo$name, function(x) {
  d <- hfile[[paste0("meta.data/", x)]]
  if ("H5D" %in% class(d)) {
    d <- d[]
  } else {
    d <- d[["values"]][]
  }
  dt <- data.table(x = d)
  setnames(dt, "x", x)
  return(dt)
})
metadata <- Reduce(cbind, metadataList)
fwrite(metadata,
       paste0("/fh/fast/gottardo_r/ytian_working/covid19_datasets/metadata/", dataset, ".tsv"),
       sep = "\t")

