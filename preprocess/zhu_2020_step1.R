library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "zhu_2020"

seu <- readRDS(file.path(data_path, dataset, "Final_nCoV_0716_upload.RDS"))
seu

### remove variable genes
VariableFeatures(seu) <- NULL

### remove the integrated assay
DefaultAssay(seu) <- "RNA"
seu[['integrated']] <- NULL

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
