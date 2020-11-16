library(SingleCellExperiment)
library(scater)
library(stringr)
library(Matrix)
library(Seurat)
library(GEOquery)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")

dataset <- "arunachalam_2020"
geo_code <- "GSE155673"

gse <- getGEO(geo_code, GSEMatrix=TRUE)
metadata <- pData(phenoData(gse[[1]]))

samples <- list.files(file.path(data_path, dataset),
                      pattern="*.mtx.gz")
samples

barcodes <- list.files(file.path(data_path, dataset),
                       pattern=".+cov\\d+_barcodes.tsv.gz")
barcodes

features <- read.delim(file.path(data_path, dataset, "GSE155673_features.tsv.gz"),
                      header = FALSE,
                      stringsAsFactors = FALSE,
                      col.names = c("ID", "Symbol", "Type"))
rownames(features) <- features$ID

### read to sce
sce_list <- vector("list", length(samples))
names(sce_list) <- str_split(samples, "_", simplify = TRUE)[, 2]

for (i in seq_along(samples)) {
    name <- names(sce_list)[[i]]
    counts <- as(readMM(file.path(data_path, dataset, samples[[i]])),
                 "dgCMatrix")
    num_cells <- dim(counts)[2]

    cell_barcodes <- read.delim(file.path(data_path, dataset, barcodes[[i]]),
                                header=FALSE)
    colnames(counts) <- cell_barcodes$V1

    coldata <- DataFrame(
        Barcode = cell_barcodes$V1,
        sample_name = rep(metadata$`sample name`[[i*2]], num_cells),
        cell_type = rep(metadata$`cell type`[[i*2]], num_cells),
        age = rep(metadata$age[[i*2]], num_cells),
        sex = rep(metadata$Sex[[i*2]], num_cells),
        disease_status = rep(metadata$disease_status[[i*2]], num_cells),
        disease_severity = rep(metadata$disease_severity[[i*2]], num_cells),
        days_since_symptom_onset = rep(metadata$days_since_symptom_onset[[i*2]],
                                       num_cells),
        row.names=NULL)

    sce <- SingleCellExperiment(list(counts = counts),
                                rowData = features,
                                colData = coldata)
    sce <- splitAltExps(sce, rowData(sce)$Type)

    counts(altExp(sce)) <- as.matrix(counts(altExp(sce)))

    rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID,
                                          rowData(sce)$Symbol)
    rownames(altExp(sce)) <- uniquifyFeatureNames(rowData(altExp(sce))$ID,
                                                  rowData(altExp(sce))$Symbol)

    sce_list[[name]] <- sce
}

all_sce <- do.call(cbind, sce_list)

seu <- CreateSeuratObject(counts(all_sce),
                          min.cells = 5)

coldata <- as.data.frame(colData(all_sce))
rownames(coldata) <- colnames(seu)
seu <- AddMetaData(seu, coldata)

seu[["ADT"]] <- CreateAssayObject(counts = counts(altExp(all_sce)))

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
