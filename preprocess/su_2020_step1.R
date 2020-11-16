library("reticulate")
library("SingleCellExperiment")
library("scater")
library("Seurat")

### import scanpy
sc <- import("scanpy")

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "su_2020"

rna <- sc$read(file.path(data_path,
                         dataset,
                         "h5ad",
                         "upto_tenth8_t1t2h_rev_P_GE_int_RAW_gex_all_samples.h5ad"))
protein <- sc$read(file.path(data_path,
                             dataset,
                             "h5ad",
                             "upto_tenth8_t1t2h_rev_P_GE_int_RAW_pro_all_samples.h5ad"))

sce <- SingleCellExperiment(
    assays      = list(counts = t(rna$X)),
    colData     = rna$obs,
    rowData     = rna$var
)
sce

sce_pro <- SingleCellExperiment(
    assays      = list(counts = t(protein$X)),
    colData     = protein$obs,
    rowData     = protein$var
)
sce_pro

altExp(sce, "ADT") <- sce_pro

seus <- vector("list", 10)

for (i in seq(1, 10)) {
    sce_batch <- sce[, sce$batch==i]
    seu <- CreateSeuratObject(counts(sce_batch),
                              meta.data=as.data.frame(colData(sce_batch)))
    seu[["ADT"]] <- CreateAssayObject(counts=counts(altExp(sce_batch)))
    seus[[i]] <- seu
}

seu <- merge(seus[[1]], unlist(seus)[2:10])
seu

saveRDS(seu, file.path(save_path,
                       paste(dataset, "step1", "rds", sep=".")))

print("DONE")

