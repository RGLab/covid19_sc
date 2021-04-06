library(Seurat)
library(SeuratDisk)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "grant_2021"

Convert(file.path(data_path, dataset, "GSE155249_supplement.h5ad"),
        dest = "h5seurat", overwrite = TRUE)
seu <- LoadH5Seurat(file.path(data_path, dataset, "GSE155249_supplement.h5seurat"))
seu

counts <- GetAssayData(seu, slot = "counts")
scale_data <- GetAssayData(seu, slot = "scale.data")
metadata <- seu@meta.data

seu2 <- CreateSeuratObject(counts = counts,
                           meta.data = metadata)
seu2 <- SetAssayData(seu2, "data", counts)
seu2 <- SetAssayData(seu2, "scale.data", scale_data)

cm <- seu2[["RNA"]]@counts
keep <- rowSums(cm>0) >= 5
print(sum(keep))
seu2 <- subset(seu2, features=rownames(seu2)[keep])
seu2

saveRDS(seu2, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
