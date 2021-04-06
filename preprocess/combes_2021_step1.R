library(SingleCellExperiment)
library(scater)
library(stringr)
library(Matrix)
library(Seurat)
library(GEOquery)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "combes_2021"
geo_code <- "GSE163668"

gse <- getGEO(geo_code, GSEMatrix=TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)

samples <- list.files(file.path(data_path, dataset),
                      pattern="*.mtx.gz")
samples

geo_accessions <- str_split(samples, "_", simplify = TRUE)[, 1]

## for (gsm in geo_accessions) {
##     files <- list.files(file.path(data_path, dataset),
##                         pattern=paste("^", gsm, sep=""),
##                         full.names=TRUE)
##     files <- files[!file.info(files)$isdir]
##     ## copy files to sub dir
##     dir.create(file.path(data_path, dataset, gsm))
##     for(f in files) {
##         oldname <- basename(f)
##         ## special case
##         if (grepl("clust1", oldname, fixed=TRUE)) {
##             oldname <- sub("\\.", "_", oldname)
##         }
##         new_name = tail(str_split(oldname, "_")[[1]], n=1)
##         file.copy(from=f,
##                   to=file.path(data_path, dataset, gsm, new_name))
##     }
## }

### read to seurat
seu_list <- vector("list", length(geo_accessions))
for (i in seq_along(geo_accessions)) {
    gsm <- geo_accessions[i]

    ## make sure matches with metadata
    print(gsm == metadata$geo_accession[[i]])

    counts <- Read10X(data.dir = file.path(data_path, dataset, gsm))
    seu <- CreateSeuratObject(counts = counts,
                              min.cells=5)
    num_cells <- dim(seu)[2]
    print(num_cells)
    sample_name <- str_split(samples[i], "_", simplify = TRUE)[, 2]
    meta <- DataFrame(
        title = rep(metadata$title[[i]], num_cells),
        cells_loaded = rep(metadata$`cells_loaded:ch1`[[i]], num_cells),
        covid_status = rep(metadata$`covid_status:ch1`[[i]], num_cells),
        phenotype = rep(metadata$`phenotype:ch1`[[i]], num_cells),
        pooled = rep(metadata$`pooled:ch1`[[i]], num_cells),
        tissue = rep(metadata$`tissue:ch1`[[i]], num_cells))

    ## pooled samples
    if (metadata$`pooled:ch1`[[i]]=="TRUE") {
        free <- read.table(file.path(data_path,
                                     dataset,
                                     gsm,
                                     "clust1.samples.txt.gz"),
                           header=TRUE, sep="\t",
                           stringsAsFactors=FALSE)
        print(dim(free))
        ## only include singlets
        sngs <- subset(free, DROPLET.TYPE=="SNG")
        seu <- seu[, sngs$BARCODE]
        meta <- meta[1:dim(seu)[2], ]
        meta$SNG.BEST.GUESS <- sngs$SNG.BEST.GUESS
        assignments <- read.table(file.path(data_path,
                                            dataset,
                                            gsm,
                                            "assignments.tsv"),
                                  header=TRUE, sep="\t")
        meta <- merge(meta,assignments, by.x="SNG.BEST.GUESS",
                      by.y="freemuxlet_cluster")
    } else {
        meta$sample <- rep(sample_name, num_cells)
    }
    rownames(meta) <- Cells(seu)
    ## seu@meta.data <- cbind(seu@meta.data, meta)
    seu <- AddMetaData(seu, metadata=as.data.frame(meta))
    seu_list[[i]] <- seu
}

##merge
seu <- merge(seu_list[[1]], y=seu_list[2:length(seu_list)])
seu

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
