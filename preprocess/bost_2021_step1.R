library(SingleCellExperiment)
library(scater)
library(stringr)
library(Matrix)
library(Seurat)
library(GEOquery)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "bost_2021"
geo_code <- "GSE157344"

gse <- getGEO(geo_code, GSEMatrix=TRUE)
metadata <- pData(phenoData(gse[[1]]))
head(metadata)

samples <- list.files(file.path(data_path, dataset),
                      pattern="*.tar.gz")
samples

gsm_to_sample <- str_split(samples, "_", simplify = TRUE)[, 2]
names(gsm_to_sample) <- str_split(samples, "_", simplify = TRUE)[, 1]
gsm_to_sample

### read to seurat
seu_list <- vector("list", length(samples))
for (i in seq_along(samples)) {
    gsm <- names(gsm_to_sample)[i]
    sample <- gsm_to_sample[[i]]
    print(sample)
    ## remove sample from file names
    folder <- file.path(data_path, dataset, "human",
                        sample, "outs", "filtered_feature_bc_matrix")

    setwd(folder)
    ## remove "_"_
    fns <- list.files(folder)
    new_fns <- str_replace(fns, "_", "")
    file.rename(fns, new_fns)
    new_fns2 <- str_replace(new_fns, paste(sample, ".", sep=""), "")
    file.rename(new_fns, new_fns2)
    counts <- Read10X(data.dir = folder)
    seu <- CreateSeuratObject(counts = counts,
                              min.cells=5)
    num_cells <- dim(seu)[2]
    print(num_cells)

    meta <- DataFrame(
        age = rep(metadata[gsm, ]$`age:ch1`, num_cells),
        clinic_status = rep(metadata[gsm, ]$`clinic status:ch1`, num_cells),
        clinical_outcome = rep(metadata[gsm, ]$`clinical outcome:ch1`, num_cells),
        sex  = rep(metadata[gsm, ]$`Sex:ch1`, num_cells),
        sofa_score = rep(metadata[gsm, ]$`sofa score:ch1`, num_cells),
        subject_id = rep(metadata[gsm, ]$`subject id:ch1`, num_cells),
        tissue = rep(metadata[gsm, ]$`tissue:ch1`, num_cells))

    rownames(meta) <- Cells(seu)
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
