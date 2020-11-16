library(Seurat)
library(stringr)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "yu_2020"

### outputs from cellranger aggr
cellranger_aggr <- file.path(data_path, dataset, "cellranger", "aggr",
                             dataset, "outs", "filtered_feature_bc_matrix")
list.files(cellranger_aggr)

data <- Read10X(data.dir = cellranger_aggr)
barcodes <- colnames(data)
barcodes[1:5]

### meatadata
aggr_csv <- read.csv(file.path(data_path, dataset, "yu_2020_aggr.csv"))
cell_orig <- str_split(barcodes, "-", simplify = TRUE)[, 2]

metadata <- aggr_csv[cell_orig, 3:6]
rownames(metadata) <- barcodes

seu = CreateSeuratObject(counts = data,
                         min.cell=5,
                         meta.data=metadata)

seu

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")

