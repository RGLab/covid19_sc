library(Seurat)

data_path <- file.path(Sys.getenv("DATA_DIR"), "data")
save_path <- file.path(Sys.getenv("DATA_DIR"), "seu")
dataset <- "schulte-schrepping_2020"

seu1 <- readRDS(file.path(data_path, dataset,
                         "seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds"))
seu1
table(seu1$orig.ident)
seu1$cohort = "cohort1"

seu2 <- readRDS(file.path(data_path, dataset,
                          "seurat_COVID19_freshWB-PBMC_cohort2_rhapsody_jonas_FG_2020-08-18.rds"))
seu2
table(seu2$orig.ident)
seu2$cohort = "cohort2"

seu <- merge(seu1, y = seu2, add.cell.ids = c("cohort1", "cohort2"),
             project = "schulte-schrepping_2020")
seu
table(seu$orig.ident)

cm <- seu[["RNA"]]@counts
keep <- rowSums(cm>0) >= 5
print(sum(keep))
seu <- subset(seu, features=rownames(seu)[keep])
seu

saveRDS(seu, file.path(save_path, paste(dataset,
                                        "step1",
                                        "rds", sep=".")))

print("DONE")
