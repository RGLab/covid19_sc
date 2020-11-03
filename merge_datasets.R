library(SeuratDisk)
library(data.table)
setDTthreads(4)

h5seuratDir <- "/fh/fast/gottardo_r/ytian_working/covid19_datasets/h5seurat"
metadataDir <- "/fh/fast/gottardo_r/ytian_working/covid19_datasets/metadata"
datasets <- c("wilk_2020",
              "chua_2020",
              "silvin_2020",
              "wen_2020",
              "lee_2020",
              "liao_2020",
              "zhu_2020",
              "arunachalam_2020",
              "yu_2020",
              "su_2020",
              "meckiff_2020")
all_meta <- fread(file.path(metadataDir, "all_meta.tsv"))

# hfile <- Connect("/fh/fast/gottardo_r/hmiller_working/silvin_2020_test.h5seurat")
# Load
seuList <- lapply(datasets, function(dataset_name) {
  message("-----", dataset_name, "-----")
  seu <- LoadH5Seurat(file.path(h5seuratDir, paste0(dataset_name, "_processed.h5seurat")),
                      assays = c("SCT", "predicted_ADT"),
                      reductions = "ref.umap",
                      meta.data = TRUE)
  # Merge appropriate meta data
  merge_field <- unique(all_meta[dataset == dataset_name]$sample_source_field)
  seu_meta <- seu@meta.data[, c(merge_field,
                                "predicted.celltype.l1",
                                "predicted.celltype.l2",
                                "predicted.celltype.l1.score",
                                "predicted.celltype.l2.score")]
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
                           "patient",
                           "sample",
                           "sex",
                           "age",
                           "disease_status",
                           "disease_severity",
                           "predicted.celltype.l1",
                           "predicted.celltype.l2",
                           "predicted.celltype.l1.score",
                           "predicted.celltype.l2.score")]
  # Set it in the correct order
  seu_meta <- seu_meta[rownames(seu@meta.data),]
  seu@meta.data <- seu_meta
  return(seu)
})

message(">>> merging datasets...")
names(seuList) <- datasets
bigSeu <- merge(seuList[[1]],
                seuList[2:length(seuList)],
                add.cell.ids = datasets,
                merge.data = TRUE,
                merge.dr = "ref.umap")
message(">>> saving to ", file.path(h5seuratDir, "covid19_datasets.h5seurat"))
SaveH5Seurat(bigSeu, file.path(h5seuratDir, "covid19_datasets.h5seurat"))
