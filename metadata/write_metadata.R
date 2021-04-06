# Extract metadata from h5seurat file and write as tsv
# Eg: Rscript write_metadata.R liao_2020

library(SeuratDisk)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

dataset <- args[1]
datadir <- file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed")
metadatadir <- file.path(Sys.getenv("DATA_DIR"), "metadata")

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

# Two samples with entries for two different sample types
if (dataset == "schulte-schrepping_2020") {
  hfile <- Connect(file.path(datadir, paste0(dataset, "_processed.h5seurat")), mode = "r+")
  sampleid <- hfile[["meta.data/sampleID"]][]
  cell_type <- hfile[["meta.data/cells"]][]
  unique_sampleid <- paste0(sampleid, "_", cell_type)
  hfile[["meta.data/sampleid_unique"]] <- unique_sampleid
  hfile$close_all()
}

# Fix Patient 17 sample which should be Patient 18.
# Add sample field with tissue to bost dataset
if (dataset == "bost_2021") {
  hfile <- Connect(file.path(datadir, paste0(dataset, "_processed.h5seurat")), mode = "r+")
  # Patient 17 with age 80 should be Patient 18.
  patient <- hfile[["meta.data/subject_id"]][]
  age <- hfile[["meta.data/age"]][]
  patient[patient == "Patient 17" & age == 80] <- "Patient 18"
  hfile[["meta.data"]]$link_delete("subject_id")
  hfile[["meta.data/subject_id"]] <- patient

  if (!"sample" %in% names(hfile[["meta.data"]])) {
    patient <- hfile[["meta.data/subject_id"]][]
    tissue <- hfile[["meta.data/tissue"]][]
    sample <- paste0(trimws(patient), "_", tissue)
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
       file.path(metadatadir, paste0(dataset, ".tsv")),
       sep = "\t")

