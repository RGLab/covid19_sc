# Merge all datasets into one giant Seurat object
# Note: R cannot hold the merged matrix in memory in sparse format
# becuase it requires indexing above R's 32-bit int limit. Thus,
# three objects are saved: bigSeu_tomerge.h5seurat, stephenson_2021_tomerge.h5seurat, covid19_datasets.h5seurat.
# covid19_datasets.h5seurat does not include the SCT assay. It is added later with merge_sct.py.

library(SeuratDisk)
library(Seurat)
library(data.table)
setDTthreads(4)

h5seuratDir <- file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed")
wholeBloodDir <- file.path(Sys.getenv("DATA_DIR"), "seu_reprocessed_singler")
metadataDir <- file.path(Sys.getenv("DATA_DIR"), "metadata")

all_meta <- fread(file.path(metadataDir, "all_meta.tsv"))
# Deal with Whole Blood
pbmc_datasets <- unique(all_meta[sample_type == "Peripheral Blood Mononuclear Cells"]$dataset)
whole_blood_datasets <- unique(all_meta[sample_type == "Whole Blood"]$dataset)

# Takes vector of gene aliases and maps to HGNC gene symbol
hgncAlias2Symbol <- fread("hgncAlias2Symbol.tsv")
mapGeneSymbol <- function(gene_aliases) {
  dt <- data.table(ALIAS = gene_aliases)
  dt <- merge(dt, hgncAlias2Symbol, all.x = TRUE, all.y = FALSE, sort = FALSE)
  dt[is.na(SYMBOL), SYMBOL := ""]
  return(dt$SYMBOL)
}

shared_genes <- readLines("genes.txt")

extract_data <- function(seu, dataset_name, assay = "SCT") {

  # Map gene alias to symbol via HUGO db (as of 12/23/2020)
  # Advice is to NOT rename features: https://github.com/satijalab/seurat/issues/1049
  # We do it only for datasets to be merged.
  seu[[assay]]@counts@Dimnames[[1]] <- mapGeneSymbol(seu[[assay]]@counts@Dimnames[[1]])
  seu[[assay]]@data@Dimnames[[1]] <- mapGeneSymbol(seu[[assay]]@data@Dimnames[[1]])
  rownames(seu[[assay]]@scale.data) <- mapGeneSymbol(rownames(seu[[assay]]@scale.data))
  seu[[assay]]@var.features <- mapGeneSymbol(seu[[assay]]@var.features)

  # Subset to only genes which successfully mapped to aliases
  seu[["SCT"]] <- subset(seu[["SCT"]], features = shared_genes)

  # Merge appropriate meta data
  merge_field <- unique(all_meta[dataset == dataset_name]$sample_source_field)
  seu_meta <- seu@meta.data[, c(merge_field,
                                "predicted.celltype.l1",
                                "predicted.celltype.l2",
                                "predicted.celltype.l3",
                                "predicted.celltype.l1.score",
                                "predicted.celltype.l2.score",
                                "predicted.celltype.l3.score",
                                "mapping.score")]
  seu_meta$id <- rownames(seu_meta)
  seu_meta$dataset <- dataset_name
  seu_meta <- merge(seu_meta,
                    all_meta,
                    by.x = c(merge_field, "dataset"),
                    by.y = c("sample", "dataset"),
                    all.x = TRUE,
                    all.y = FALSE)
  seu_meta$sample <- seu_meta[[merge_field]]
  rownames(seu_meta) <- seu_meta$id
  seu_meta <- seu_meta[, c("dataset",
                           "tissue",
                           "sample_type",
                           "sample_type_note",
                           "patient",
                           "sample",
                           "race_reported",
                           "sex",
                           "age",
                           "disease_status",
                           "disease_severity",
                           "days_since_symptom_onset",
                           "predicted.celltype.l1",
                           "predicted.celltype.l2",
                           "predicted.celltype.l3",
                           "predicted.celltype.l1.score",
                           "predicted.celltype.l2.score",
                           "predicted.celltype.l3.score",
                           "mapping.score")]
  # Set it in the correct order
  seu_meta <- seu_meta[rownames(seu@meta.data),]
  seu@meta.data <- seu_meta
  DefaultAssay(seu) <- "SCT"
  seu[["ref.umap"]]@assay.used <- "SCT"
  DietSeurat(seu,
             assays = c("SCT", "predicted_ADT"),
             dimreducs = c("ref.umap", "ref.spca"))
}


# Create list of all of them except stephenson. Save stephenson to a separate file.

# Just PBMCs
seuList_pbmc <- lapply(pbmc_datasets, function(dataset_name) {
  message("-----", dataset_name, "-----")
  loadDir <- h5seuratDir

  seu <- LoadH5Seurat(file.path(loadDir, paste0(dataset_name, "_processed.h5seurat")),
                      assays = c("SCT", "predicted_ADT"),
                      reductions = c("ref.umap", "ref.spca"),
                      meta.data = TRUE,
                      verbose = FALSE)

  # Subset to only PBMCs
  pbmc_samples <- all_meta[dataset == dataset_name &
                             sample_type == "Peripheral Blood Mononuclear Cells"]$sample
  sample_field <- unique(all_meta[dataset == dataset_name]$sample_source_field)
  seu <- seu[,seu[[sample_field]][,] %in% pbmc_samples]

  return(extract_data(seu, dataset_name, "SCT"))
})

names(seuList_pbmc) <- pbmc_datasets


# Just whole blood
seuList_wb <- lapply(whole_blood_datasets, function(dataset_name) {
  message("----- ", dataset_name, " whole blood -----")
  path <- file.path(wholeBloodDir, list.files(wholeBloodDir, pattern = dataset_name))
  if (length(path) != 1) stop("Found ", length(path), " whole blood files for ", dataset_name)
  seu <- readRDS(path)
  seu <- seu[, !seu$monaco_main %in% c("Neutrophils", "Basophils")]

  # For 2 datasets where it is missing, add derived "sample" column
  # (normally would happen when writing metadata)
  if (dataset_name == "silvin_2020") {
    seu <- AddMetaData(seu,
                       paste0(trimws(seu$Characteristics.individual.),
                              "_",
                              seu$Factor.Value.sampling.time.point.),
                       "sample")
  } else if (dataset_name == "bost_2021") {
    seu <- AddMetaData(seu,
                       paste0(trimws(seu$subject_id),
                              "_",
                              seu$tissue),
                       "sample")
  } else if (dataset_name == "schulte-schrepping_2020") {
    seu <- AddMetaData(seu,
                       paste0(seu$sampleID, "_", seu$cells),
                       "sampleid_unique")
  }

  return(extract_data(seu, dataset_name, "SCT"))
})

names(seuList_wb) <- whole_blood_datasets

# Now merge all but stephenson, with SCT. save it.
message(">>> saving stephenson_2021_tomerge to ", h5seuratDir)
stephenson_tomerge <- seuList_pbmc[["stephenson_2021"]]
stephenson_tomerge <- RenameCells(stephenson_tomerge, add.cell.id = "stephenson_2021")
SaveH5Seurat(object = stephenson_tomerge,
             filename = file.path(h5seuratDir, "stephenson_2021_tomerge.h5seurat"),
             verbose = TRUE,
             overwrite = TRUE)
seuList_pbmc[["stephenson_2021"]] <- NULL
seuList <- c(seuList_pbmc, seuList_wb)

message(">>> Merging all but stephenson...")
bigSeu <- merge(seuList[[1]],
                seuList[2:length(seuList)],
                add.cell.ids = names(seuList),
                merge.data = TRUE,
                merge.dr = c("ref.umap", "ref.spca"))

message(">>> saving bigSeu...")
SaveH5Seurat(bigSeu,
             file.path(h5seuratDir, "bigSeu_tomerge.h5seurat"),
             verbose = TRUE,
             overwrite = TRUE)

message(">>> Merging stephenson ADT...")
# Merge in stephenson, excluding SCT.
DefaultAssay(bigSeu) <- "predicted_ADT"
bigSeu[["ref.spca"]]@global <- TRUE
bigSeu <- DietSeurat(bigSeu,
                     assays = "predicted_ADT",
                     dimreducs = c("ref.umap", "ref.spca"))
DefaultAssay(stephenson_tomerge) <- "predicted_ADT"
stephenson_tomerge[["ref.spca"]]@global <- TRUE
stephenson_tomerge <- DietSeurat(stephenson_tomerge,
                             assays = "predicted_ADT",
                             dimreducs = c("ref.umap", "ref.spca"))
bigSeu <- merge(bigSeu,
                stephenson_tomerge,
                merge.data = TRUE,
                merge.dr = c("ref.umap", "ref.spca"))

bigSeu[["ref.umap"]]@assay.used <- "SCT"
message(">>> saving to ", file.path(h5seuratDir, "covid19_datasets.h5seurat"))
SaveH5Seurat(bigSeu, file.path(h5seuratDir, "covid19_datasets.h5seurat"), overwrite = TRUE)
# Use Python to manually merge SCT and add it to the file.

