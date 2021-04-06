library(Matrix)
library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "heming_2021"

data <- Read10X(file.path(data_path, dataset))
seu = CreateSeuratObject(counts = data,
                         min.cells=5)
seu

cluster <- read.csv(file.path(data_path, dataset, "GSE163005_annotation_cluster.csv"),
                    col.names=c("barcode", "celltype"))
head(cluster)

dx  <- read.csv(file.path(data_path, dataset, "GSE163005_annotation_dx.csv"),
                    col.names=c("barcode", "diagnosis"))
head(dx)

t_cluster <- read.csv(file.path(data_path, dataset, "GSE163005_annotation_tcells_cluster.csv"),
                    col.names=c("barcode", "tcells_cluster"))
head(t_cluster)

dim(cluster)
dim(dx)
dim(t_cluster)

anno <- merge(cluster, dx, by="barcode", all=TRUE)
anno <- merge(anno, t_cluster, by="barcode", all=TRUE)
rownames(anno) <- anno$barcode
head(anno)
dim(anno)

## subset
seu <- seu[, anno$barcode]
seu <- AddMetaData(seu, anno)
colnames(seu[[]])

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
