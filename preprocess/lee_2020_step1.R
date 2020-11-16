library(SingleCellExperiment)
library(scater)
library(stringr)
library(Matrix)
library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "lee_2020"

features <- read.delim(file.path(data_path,
                                 dataset,
                                 "GSE149689_features.tsv.gz"),
                      header = FALSE,
                      stringsAsFactors = FALSE,
                      col.names = c("ID", "Symbol", "Type"))
rownames(features) <- features$ID
head(features)

counts <- as(readMM(file.path(data_path, dataset, "GSE149689_matrix.mtx.gz")),
             "dgCMatrix")

cell_barcodes <- read.delim(file.path(data_path,
                                      dataset,
                                      "GSE149689_barcodes.tsv.gz"),
                            header=FALSE)
head(cell_barcodes)

"The suffix of each barcode sequence represents patient information.
'-1' - 'Sample1' - 'nCoV 1 scRNA-seq'
'-2' - 'Sample2' - 'nCoV 2 scRNA-seq'
'-3' - 'Sample3' - 'Flu 1 scRNA-seq'
'-4' - 'Sample4' - 'Flu 2 scRNA-seq'
'-5' - 'Sample5' - 'Normal 1 scRNA-seq'
'-6' - 'Sample6' - 'Flu 3 scRNA-seq'
'-7' - 'Sample7' - 'Flu 4 scRNA-seq'
'-8' - 'Sample8' - 'Flu 5 scRNA-seq'
'-9' - 'Sample9' - 'nCoV 3 scRNA-seq'
'-10' - 'Sample10' - 'nCoV 4 scRNA-seq'
'-11' - 'Sample11' - 'nCoV 5 scRNA-seq'
'-12' - 'Sample12' - 'nCoV 6 scRNA-seq'
'-13' - 'Sample13' - 'Normal 2 scRNA-seq'
'-14' - 'Sample14' - 'Normal 3 scRNA-seq'
'-15' - 'Sample15' - 'nCoV 7 scRNA-seq'
'-16' - 'Sample16' - 'nCoV 8 scRNA-seq'
'-17' - 'Sample17' - 'nCoV 9 scRNA-seq'
'-18' - 'Sample18' - 'nCoV 10 scRNA-seq'
'-19' - 'Sample19' - 'Normal 4 scRNA-seq'
'-20' - 'Sample20' - 'nCoV 11 scRNA-seq'"

coldata <- DataFrame(Barcode = cell_barcodes$V1,
                     row.names=cell_barcodes$V1)
colnames(counts) <- cell_barcodes$V1

samples <- paste0("Sample", seq(1, 20), sep="")
samples[1]

patients <- c('nCoV 1 scRNA-seq',
              'nCoV 2 scRNA-seq',
              'Flu 1 scRNA-seq',
              'Flu 2 scRNA-seq',
              'Normal 1 scRNA-seq',
              'Flu 3 scRNA-seq',
              'Flu 4 scRNA-seq',
              'Flu 5 scRNA-seq',
              'nCoV 3 scRNA-seq',
              'nCoV 4 scRNA-seq',
              'nCoV 5 scRNA-seq',
              'nCoV 6 scRNA-seq',
              'Normal 2 scRNA-seq',
              'Normal 3 scRNA-seq',
              'nCoV 7 scRNA-seq',
              'nCoV 8 scRNA-seq',
              'nCoV 9 scRNA-seq',
              'nCoV 10 scRNA-seq',
              'Normal 4 scRNA-seq',
              'nCoV 11 scRNA-seq')

keys <- coldata$Barcode %>%
    strsplit("-") %>%
    sapply("[", 2 ) %>%
    as.integer()

coldata$Sample <- samples[keys]
coldata$Patient <- patients[keys]
head(coldata)

sce <- SingleCellExperiment(list(counts = counts),
                            rowData = features,
                            colData = coldata)

rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID,
                                      rowData(sce)$Symbol)
sce

seu <- CreateSeuratObject(counts(sce),
                          min.cells = 5)

coldata <- as.data.frame(colData(sce))
seu <- AddMetaData(seu, coldata)

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
