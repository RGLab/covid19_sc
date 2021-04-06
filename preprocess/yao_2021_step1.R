library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "yao_2021"

samples <- list.files(file.path(data_path, dataset),
                      pattern="*.h5")
samples

## reorder samples
reorder_index <- c(1, 2, 3, 7, 8, 9, 4, 5, 6)
samples <- samples[order(reorder_index)]
samples

seu_list <- vector("list", length(samples))
for (i in seq_along(samples)) {
    data <- Read10X_h5(file.path(data_path, dataset, samples[[i]]))
    current <- CreateSeuratObject(data)
    current <- RenameCells(current,
                           new.names = paste(colnames(x = current[["RNA"]]), i, sep="_"))
    seu_list[[i]] <- current
}
length(seu_list)

seu <- merge(seu_list[[1]], y = seu_list[2:9])
seu <- CreateSeuratObject(GetAssayData(seu, slot="counts"),
                          min.cells = 5,
                          project="yao_2021")
seu

anno <- read.csv(file.path(data_path, dataset, "GSE154567_annotation.csv"),
                 row.names=1)
dim(anno)

all(rownames(anno) %in% colnames(seu))

## subsetting
seu <- seu[, rownames(anno)]
seu

seu <- AddMetaData(seu, anno)

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
