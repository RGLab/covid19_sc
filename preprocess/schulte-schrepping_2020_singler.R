library(Seurat)
library(SingleR)

data_path <- Sys.getenv("DATA_DIR")
dataset <- "schulte-schrepping_2020"

seu <- readRDS(file.path(data_path, "seu_reprocessed",
                         paste(dataset, "processed.rds", sep="_")))
seu

##subset wholeblood
table(seu$cells)

seu <- subset(seu, subset = cells == "Whole_blood")
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
