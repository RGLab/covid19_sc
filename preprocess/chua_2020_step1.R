library(Seurat)

### The _main file contains all samples from the nasopharynx, while the _loc file contains data from nasopharyngeal, protected specimen brush, and bronchial lavage samples of two patients.

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "chua_2020"

seu_loc <- readRDS(file.path(data_path, dataset, "covid_nbt_loc.rds"))
seu_main <- readRDS(file.path(data_path, dataset, "covid_nbt_main.rds"))

seu_loc
seu_main

### remove ns samples from loc, which are already in main
'%!in%' <- function(x,y)!('%in%'(x,y))
seu_loc <- subset(seu_loc, subset = sample %!in% c("BIH-CoV-01_NS_1",
                                                  "BIH-CoV-04_NS_1"))
seu_loc

### merge
seu <- merge(seu_loc, y = seu_main,
             add.cell.ids = c("loc", "main"))
seu

### filter genes
cm <- seu[["RNA"]]@counts
keep <- rowSums(cm>0) >= 5
seu <- subset(seu, features=keep)
seu

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
