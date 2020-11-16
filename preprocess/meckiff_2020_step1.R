library(Matrix)
library(Seurat)

### Read Counts
data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "meckiff_2020"

cd4t0n6_coldata <- read.delim(file.path(data_path,
                                        dataset,
                                        "GSE152522_cd4t0n6_annotation.txt.gz"))
cd4t0n6_counts <- read.delim(file.path(data_path,
                                       dataset,
                                       "GSE152522_cd4t0n6_umi.txt.gz"),
                             row.names=1)

dim(cd4t0n6_coldata)
dim(cd4t0n6_counts)

cd4t24_coldata <- read.delim(file.path(data_path,
                                        dataset,
                                        "GSE152522_cd4t24_annotation.txt.gz"))
cd4t24_counts <- read.delim(file.path(data_path,
                                      dataset,
                                      "GSE152522_cd4t24_umi.txt.gz"),
                            row.names=1)

dim(cd4t24_coldata)
dim(cd4t24_counts)

cd4t6_coldata <- read.delim(file.path(data_path,
                                        dataset,
                                        "GSE152522_cd4t6_annotation.txt.gz"))
cd4t6_counts <- read.delim(file.path(data_path,
                                     dataset,
                                     "GSE152522_cd4t6_umi.txt.gz"),
                           row.names=1)

dim(cd4t6_coldata)
dim(cd4t6_counts)

cd4t0n6_counts <- as(as.matrix(cd4t0n6_counts), "dgCMatrix")
cd4t24_counts <- as(as.matrix(cd4t24_counts), "dgCMatrix")
cd4t6_counts <- as(as.matrix(cd4t6_counts), "dgCMatrix")

rownames(cd4t0n6_coldata) <- colnames(cd4t0n6_counts)
rownames(cd4t24_coldata) <- colnames(cd4t24_counts)
rownames(cd4t6_coldata) <- colnames(cd4t6_counts)

seu_t0n6 <- CreateSeuratObject(counts = cd4t0n6_counts,
                               meta.data = cd4t0n6_coldata)
seu_t24 <- CreateSeuratObject(counts = cd4t24_counts,
                              meta.data = cd4t24_coldata)
seu_t6 <- CreateSeuratObject(counts = cd4t6_counts,
                             meta.data = cd4t6_coldata)

### merge
seu <- merge(seu_t0n6, y = seu_t24,
             add.cell.ids = c("t0n6", "t24"))
seu

### remove hashtag nas
seu <- subset(seu, subset = orig.HT_ID.global=="Singlet")
seu

### filter genes
cm <- seu[["RNA"]]@counts
keep <- rowSums(cm>0) >= 5
seu <- subset(seu, features=keep)
seu

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))
print("DONE")
