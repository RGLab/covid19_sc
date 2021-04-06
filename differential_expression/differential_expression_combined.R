# Get differentially expressed marker genes for L1, L2, and L3 cell types for the combined dataset

message(getwd())
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(HDF5Array)
library(MAST)
library(slurmR)

debugDir <- file.path(Sys.getenv("DEBUG_DIR"), "covid_DE", "combined")
h5seuratDir <- file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed")
h5SCEDir <- file.path(h5seuratDir, "combined_sce")

message(">>> debugDir: ", debugDir)
message(">>> h5seuratDir: ", h5seuratDir)
message(">>>h5SCEDir: ", h5SCEDir)

# h5SCEDir <- "/fh/scratch/delete30/gottardo_r/hmiller/example_datasets/combined_sce"
# h5seuratDir <- "/fh/scratch/delete30/gottardo_r/hmiller/example_datasets"

message(">>> Getting cell types from ", file.path(h5seuratDir, "covid19_datasets.h5seurat"))
seu <- Connect(file.path(h5seuratDir, "covid19_datasets.h5seurat"))
celltypesL1 <- table(seu[["meta.data"]][["predicted.celltype.l1"]][])
celltypesL1 <- names(celltypesL1)[celltypesL1 >5]
celltypesL2 <- table(seu[["meta.data"]][["predicted.celltype.l2"]][])
celltypesL2 <- names(celltypesL2)[celltypesL2 > 5]
celltypesL3 <- table(seu[["meta.data"]][["predicted.celltype.l3"]][])
celltypesL3 <- names(celltypesL3)[celltypesL3 > 5]
seu$close_all()

slurmR_tmp_dir <- file.path(Sys.getenv("DEBUG_DIR"), "covid_DE", "slurmr")
if (!dir.exists(slurmR_tmp_dir)) dir.create(slurmR_tmp_dir)

message(">>> slurmR_tmp_dir: ", slurmR_tmp_dir)

FindMarkersForCelltype <- function(celltype,
                                   inpath,
                                   label) {
  sce <- loadHDF5SummarizedExperiment(inpath, prefix = "combined_")

  cells.1 <- Cells(sce)[which(sce[[label]] == celltype)]
  cells.2 <- Cells(sce)[which(sce[[label]] != celltype)]
  fc.results <- Seurat::FoldChange(logcounts(sce),
                           cells.1 = cells.1,
                           cells.2 = cells.2,
                           fc.name = "avg_log2FC",
                           mean.fxn = function(x) {
                             return(log(x = rowMeans(x = expm1(x = x)) + 1, base = 2))
                           })
  markers <- Seurat::FindMarkers(logcounts(sce),
                           cells.1 = cells.1,
                           cells.2 = cells.2,
                           only.pos = TRUE,
                           min.pct = 0.25,
                           logfc.threshold = 0.25,
                           test.use = "MAST",
                           fc.results = fc.results)
  markers$gene <- rownames(markers)
  markers$cluster <- celltype
  return(markers)
}

l1_job <- Slurm_lapply(celltypesL1, FindMarkersForCelltype,
                       inpath = h5SCEDir,
                       label = "predicted.celltype.l1",

                       njobs = length(celltypesL1),
                       mc.cores = 1,
                       job_name = paste0("DE_l1_combined_", as.character(Sys.Date())),
                       tmp_path = slurmR_tmp_dir,
                       sbatch_opt = list("constraint" = "gizmok",
                                         "cpus-per-task" = "16"),
                       seeds = 1:length(celltypesL1),
                       plan = "submit",
                       overwrite = TRUE)


l2_job <- Slurm_lapply(celltypesL2, FindMarkersForCelltype,
                       inpath = h5SCEDir,
                       label = "predicted.celltype.l2",

                       njobs = length(celltypesL2),
                       mc.cores = 1,
                       job_name = paste0("DE_l2_combined_", as.character(Sys.Date())),
                       tmp_path = slurmR_tmp_dir,
                       sbatch_opt = list("constraint" = "gizmok",
                                         "cpus-per-task" = "16"),
                       seeds = 1:length(celltypesL2),
                       plan = "submit",
                       overwrite = TRUE)

l3_job <- Slurm_lapply(celltypesL3, FindMarkersForCelltype,
                       inpath = h5SCEDir,
                       label = "predicted.celltype.l3",

                       njobs = length(celltypesL3),
                       mc.cores = 1,
                       job_name = paste0("DE_l3_combined_", as.character(Sys.Date())),
                       tmp_path = slurmR_tmp_dir,
                       sbatch_opt = list("constraint" = "gizmok",
                                         "cpus-per-task" = "16"),
                       seeds = 1:length(celltypesL3),
                       plan = "submit",
                       overwrite = TRUE)

collectResults <- function(de_job) {
  still_running <- ( length(status(de_job)$running) + length(status(de_job)$pending) )> 0
  while(still_running) {
    Sys.sleep(1)
    still_running <- (length(status(de_job)$running) + length(status(de_job)$pending)) > 0
  }

  markerList <- Slurm_collect(de_job)
  markerdf <- rbindlist(markerList)
  setDT(markerdf)
  setorder(markerdf, p_val_adj, -avg_log2FC)

  # Return the top 30 markers with p-value < 0.05
  markerdf[p_val_adj < 0.05, head(.SD, 30), cluster]
}

Sys.sleep(5)

message(">>> Waiting for results for l1...")
de_markers_l1 <- collectResults(l1_job)
message(">>> Collected results for l1")
message(">>> Writing l1 results to ", file.path(debugDir, "de_markers_l1.csv"))
fwrite(de_markers_l1, file.path(debugDir, "de_markers_l1.csv"))

message(">>> Waiting for results for l2...")
de_markers_l2 <- collectResults(l2_job)
message(">>> Collected results for l2")
message(">>> Writing l2 results to ", file.path(debugDir, "de_markers_l2.csv"))
write.csv(de_markers_l2, file.path(debugDir, "de_markers_l2.csv"))

message(">>> Waiting for results for l3...")
de_markers_l3 <- collectResults(l3_job)
message(">>> Collected results for l3")
message(">>> Writing l3 results to ", file.path(debugDir, "de_markers_l3.csv"))
write.csv(de_markers_l3, file.path(debugDir, "de_markers_l3.csv"))


# Save results to hdf5 file in @misc slot
h5seu <- Connect(file.path(h5seuratDir, "covid19_datasets.h5seurat"), mode = "r+")

h5seu[["misc"]]$create_group("de_results")
h5seu[["misc"]][["de_results"]][["de_markers_L1"]] <- de_markers_l1
h5seu[["misc"]][["de_results"]][["de_markers_L2"]] <- de_markers_l2
h5seu[["misc"]][["de_results"]][["de_markers_L3"]] <- de_markers_l3

h5seu$close_all()
