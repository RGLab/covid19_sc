library(Matrix)
library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "kusnadi_2021"

coldata <- read.delim(file.path(data_path,
                                dataset,
                                "GSE153931_cd8_t24_processed_data_annotations.txt.gz"))

counts <- read.csv(file.path(data_path,
                               dataset,
                               "GSE153931_cd8_t24_processed_data_umi_counts.txt.gz"),
                   row.names=1)

dim(coldata)
dim(counts)

rownames(coldata) <- colnames(counts)

seu <- CreateSeuratObject(counts = counts,
                          meta.data = coldata,
                          min.cells=5)
seu

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))
print("DONE")
