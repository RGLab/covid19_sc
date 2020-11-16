library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "wilk_2020"

seu <- readRDS(file.path(data_path, dataset, "blish_covid.seu.rds"))
seu

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
