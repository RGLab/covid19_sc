library(Seurat)
library(SingleR)

data_path <- "/fh/fast/gottardo_r/ytian_working/covid19_datasets"
dataset <- "bost_2021"

seu <- readRDS(file.path(data_path, "seu_reprocessed",
                         paste(dataset, "processed.rds", sep="_")))
seu

##subset wholeblood
table(seu$tissue)

seu <- subset(seu, subset = tissue == "Blood")
seu

ref <- MonacoImmuneData()
monaco_main <- SingleR(test=as.SingleCellExperiment(seu),
                       ref=ref,
                       labels=ref$label.main)
seu$monaco_main <- factor(monaco_main$pruned.labels)

table(seu$monaco_main)

meta <- seu@meta.data
write.csv(meta, file.path(data_path, "csv", paste(dataset,
                                                  "wholeblood",
                                                  "metadata",
                                                  "csv",
                                                  sep=".")))
          
saveRDS(seu, file.path(data_path, "seu_reprocessed_singler",
                       paste(dataset, "wholeblood_processed.rds", sep="_")))

print("DONE")
