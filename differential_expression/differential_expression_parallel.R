# Performs differential expression analysis for each of the L3 cell types
# based on Seurat PBMC reference using the MAST test.
# Results will be saved in misc slot in h5seurat file.
#
# Uses slurmR to launch separate slurm jobs for each cluster.
#
# ----- Parse options ------
message(getwd())
library(optparse)
option_list = list(
  make_option(c("-D", "--dataset"),
              help = "name of dataset",
              default = "silvin_2020"),
  make_option(c("-i", "--indir"),
              help = "directory with preprocessed h5seurat datasets",
              default = file.path(Sys.getenv("DATA_DIR"), "h5seurat_reprocessed")),
  make_option(c("-d", "--debugdir"),
              help = "directory to store final processed dataset",
              default = file.path(Sys.getenv("DEBUG_DIR"), "covid_DE")),
  make_option(c("-c", "--corespertask"),
              help = "number of cores to allocate to each sbatch job",
              default = "4")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
message(paste0(capture.output(opt), collaps = "\n"))
# check args
inpath <- file.path(opt$indir, paste0(opt$dataset, "_processed.h5seurat"))
if (!file.exists(inpath)) {
  stop("could not find ", inpath)
}
debugDir <- file.path(opt$debugdir, opt$dataset)
if (!dir.exists(opt$debugdir)) {
  stop("could not find ", opt$debugdir)
}
if (!dir.exists(debugDir)) dir.create(debugDir)
dataset <- opt$dataset
# ----- START -----
# https://satijalab.org/seurat/articles/de_vignette.html

library(Seurat)
library(SeuratDisk)
library(slurmR)
library(data.table)

seu <- Connect(inpath)
celltypesL1 <- table(seu[["meta.data"]][["predicted.celltype.l1"]][])
celltypesL1 <- names(celltypesL1)[celltypesL1 >5]
celltypesL2 <- table(seu[["meta.data"]][["predicted.celltype.l2"]][])
celltypesL2 <- names(celltypesL2)[celltypesL2 > 5]
celltypesL3 <- table(seu[["meta.data"]][["predicted.celltype.l3"]][])
celltypesL3 <- names(celltypesL3)[celltypesL3 > 5]
assay <- ifelse("SCT" %in% names(seu[["assays"]]), "SCT", "RNA")
seu$close_all()

slurmR_tmp_dir <- file.path(opt$debugdir, "slurmr")
if (!dir.exists(slurmR_tmp_dir)) dir.create(slurmR_tmp_dir)

FindMarkersForCelltype <- function(celltype,
                                   inpath,
                                   label,
                                   assay) {
  seu <- LoadH5Seurat(inpath,
                      assays = "SCT")
  Idents(seu) <- label

  markers <- FindMarkers(object = seu,
              assay = assay,
              slot = "data",
              ident.1 = celltype,
              only.pos = TRUE,
              min.pct = 0.25,
              logfc.threshold = 0.25,
              test.use = "MAST")
  markers$gene <- rownames(markers)
  markers$cluster <- celltype
  return(markers)
}


l1_job <- Slurm_lapply(celltypesL1, FindMarkersForCelltype,
                       inpath = inpath,
                       label = "predicted.celltype.l1",
                       assay = assay,

                       njobs = length(celltypesL1),
                       mc.cores = 1,
                       job_name = paste0("DE_l1_", opt$dataset, "_", as.character(Sys.Date())),
                       tmp_path = slurmR_tmp_dir,
                       sbatch_opt = list("constraint" = "gizmok",
                                         "cpus-per-task" = opt$corespertask),
                       seeds = 1:length(celltypesL1),
                       plan = "submit",
                       overwrite = TRUE)


l2_job <- Slurm_lapply(celltypesL2, FindMarkersForCelltype,
                       inpath = inpath,
                       label = "predicted.celltype.l2",
                       assay = assay,

                       njobs = length(celltypesL2),
                       mc.cores = 1,
                       job_name = paste0("DE_l2_", opt$dataset, "_", as.character(Sys.Date())),
                       tmp_path = slurmR_tmp_dir,
                       sbatch_opt = list("constraint" = "gizmok",
                                         "cpus-per-task" = opt$corespertask),
                       seeds = 1:length(celltypesL2),
                       plan = "submit",
                       overwrite = TRUE)

l3_job <- Slurm_lapply(celltypesL3, FindMarkersForCelltype,
                       inpath = inpath,
                       label = "predicted.celltype.l3",
                       assay = assay,

                       njobs = length(celltypesL3),
                       mc.cores = 1,
                       job_name = paste0("DE_l3_", opt$dataset, "_", as.character(Sys.Date())),
                       tmp_path = slurmR_tmp_dir,
                       sbatch_opt = list("constraint" = "gizmok",
                                         "cpus-per-task" = opt$corespertask),
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
h5seu <- Connect(inpath, mode = "r+")

h5seu[["misc"]]$create_group("de_results")
h5seu[["misc"]][["de_results"]][["de_markers_L1"]] <- de_markers_l1
h5seu[["misc"]][["de_results"]][["de_markers_L2"]] <- de_markers_l2
h5seu[["misc"]][["de_results"]][["de_markers_L3"]] <- de_markers_l3

h5seu$close_all()
