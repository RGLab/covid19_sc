library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "liao_2020"

seu <- readRDS(file.path(data_path, dataset, "nCoV.rds"))
seu

### remove variable genes
VariableFeatures(seu) <- NULL

### remove the integrated assay
seu[['integrated']] <- NULL

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
