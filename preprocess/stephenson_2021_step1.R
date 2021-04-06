library(Seurat)
library(SingleCellExperiment)
library(reticulate)
library(Matrix)

### import scanpy
sc <- import("scanpy")

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "stephenson_2021"

rna <- sc$read(file.path(data_path,
                         dataset,
                         "rna.h5ad"))
rna

adt <- sc$read(file.path(data_path,
                         dataset,
                         "adt.h5ad"))
adt

sce <- SingleCellExperiment(
    assays      = list(counts = t(rna$X)),
    colData     = rna$obs,
    rowData     = rna$var
)
sce

sce_adt <- SingleCellExperiment(
    assays      = list(counts = t(adt$X)),
    colData     = adt$obs,
    rowData     = adt$var
)
sce_adt

altExp(sce, "ADT") <- sce_adt

seus <- vector("list", 3)
sites <- levels(sce$Site)
sites

for (i in seq(1, 3)) {
    sce_batch <- sce[, sce$Site==sites[i]]
    seu <- CreateSeuratObject(counts(sce_batch),
                              meta.data=as.data.frame(colData(sce_batch)))
    seu[["ADT"]] <- CreateAssayObject(counts=counts(altExp(sce_batch)))
    seus[[i]] <- seu
}

seu <- merge(seus[[1]], seus[2:3])
seu

### filter genes
cm <- seu[["RNA"]]@counts
keep <- rowSums(cm>0) >= 5
print(sum(keep))
seu2 <- subset(seu, features=rownames(seu)[keep])
seu2[["ADT"]] <- seu[["ADT"]]
seu2

saveRDS(seu2, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
