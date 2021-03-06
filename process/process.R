# Process seurat dataset. Save output as rds and h5seurat object.
# Intermediate objects will be saved to debugdir
#
#
# ----- Parse options ------
message(getwd())
library(optparse)
option_list = list(
  make_option(c("-D", "--dataset"),
              help = "name of dataset"),
  make_option(c("-i", "--indir"),
              help = "directory with preprocessed datasets",
              default = file.path(Sys.getenv("DATA_DIR"), "seu")),
  make_option(c("-o", "--outdir"),
              help = "directory to store final processed dataset",
              default = file.path(Sys.getenv("DATA_DIR"), "seu")),
  make_option(c("-d", "--debugdir"),
              help = "path to store intermediate objects for debugging",
              default = Sys.getenv("DEBUG_DIR")),
  make_option(c("-r", "--referencedat"),
              help = "path to reference data for mapping",
              default = Sys.getenv("REFERENCE_DATASET")),
  make_option(c("-s", "--h5dir"),
              help = "directory to save h5Seurat object",
              default = file.path(Sys.getenv("DATA_DIR", "h5seurat"))),
  make_option(c("-l", "--plotdir"),
              help = "path to write plots",
              default = file.path(Sys.getenv("DATA_DIR", "plots_seurat"))),
  make_option(c("-b", "--batch"),
              help = "(optional) batch field. If specified, batches will be noramalized separately and aligned using SCTransform methods.",
              default = "null"),
  make_option(c("-n", "--normalized"),
              help = "Is the dataset already normalized? If specified, skip normalization step",
              action = "store_true",
              default = FALSE)
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
message(paste0(capture.output(opt), collaps = "\n"))
# check args
inpath <- file.path(opt$indir, paste0(opt$dataset, ".step1.rds"))
if (!file.exists(inpath)) {
  stop("could not find ", inpath)
}
outpath_seu <- file.path(opt$outdir, paste0(opt$dataset, "_processed.rds"))
if (!dir.exists(opt$outdir)) {
  stop("could not find ", opt$outdir)
}
outpath_h5s <- file.path(opt$h5dir, paste0(opt$dataset, "_processed.h5seurat"))
if (!dir.exists(opt$h5dir)) {
  stop("could not find ", opt$h5dir)
}
debugdir <- file.path(opt$debugdir, opt$dataset)
if (!dir.exists(debugdir)) {
  message("creating debug dir: ", debugdir)
  dir.create(debugdir)
}
plotdir <- file.path(opt$plotdir, opt$dataset)
if (!dir.exists(plotdir)) {
  dir.create(plotdir)
}

if (!file.exists(opt$referencedat)) {
  stop("could not find ", opt$referencedat)
}

# ----- Setup -----

# Load functions and libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(magrittr)
  library(SeuratDisk)
  library(ggplot2)
})

set.seed(1001)
ncores <- Sys.getenv("SLURM_CPUS_PER_TASK")
# Set the maxSize for future to 2G per core
options(future.globals.maxSize = 2 * as.integer(ncores) * 1024^3)
if (FALSE) {
  # if (ncores != "") {
  # TODO: Set reproducible seed
  # https://github.com/satijalab/seurat/issues/3622
  message(">>> Setting future strategy with ", ncores, " cores")
  future::plan(strategy = "multicore", workers = as.integer(ncores))
}

# ----- START -----
message(">>> Reading input: ", inpath)
seu <- readRDS(inpath)
seu <- UpdateSeuratObject(seu)
DefaultAssay(seu) <- "RNA"



### ----- Process RNA -----
message(">>> Processing RNA...")
if (opt$normalized) {

  # For studies wih normalized RNA data but not counts,
  # calculate variable features, PCA, clusters, UMAP based on
  # provided normalized values first, labeling these results
  # as "pub" to indicated that is based on the published
  # normalized values. We also run sctransform on these datasets
  # so that they are comparable to other datasets.

  message(">>> No counts data. Skipping QC...")
  message(">>> Running standard analyses on normalized RNA data...")
  if (length(VariableFeatures(seu)) == 0) {
    message(">>> Finding variable features...")
    seu <- FindVariableFeatures(seu,
                                nfeatures = 3000)
  }
  if (sum(dim(seu[[DefaultAssay(seu)]]@scale.data)) == 0) {
    message(">>> Scaling data")
    seu <- ScaleData(seu)
  }

  message(">>> Calculating PCs...")
  seu <- RunPCA(seu,
                reduction.name = "pca.pub",
                verbose=TRUE)

  ### clustering based on RNA
  message(">>> finding clusters...")
  seu <- FindNeighbors(seu,
                       reduction = "pca.pub",
                       dims = 1:30,
                       graph.name="rna.snn.pub") %>%
    FindClusters(graph.name="rna.snn.pub")
  seu[["rnaClusterID.pub"]] <- Idents(seu)

  ### RNA UMAP
  message(">>> Calculating RNA UMAP")
  seu <- RunUMAP(seu,
                 reduction = 'pca.pub',
                 dims = 1:30,
                 assay = 'RNA',
                 reduction.name = 'rna.umap.pub',
                 reduction.key = 'rnaUMAP_')

} else {

  message(">>> QC cells...")
  # Objects which were converted from non-Seurat object may not have nCount_RNA and nFeature_RNA.
  if (!"nCount_RNA" %in% names(seu@meta.data)) {
    nCount = colSums(seu, slot = "counts")  # nCount_RNA
    nFeature = colSums(x = GetAssayData(seu, slot = "counts") > 0)  # nFeatureRNA
    seu <- AddMetaData(seu, data.frame(nCount_RNA = nCount, nFeature_RNA = nFeature))
  }
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  # Write out some diagnostic plots
  pdf(file=file.path(plotdir, "QC.pdf"),
      paper="USr")
  plot1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept = 50)
  plot3 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept = 200)
  plot1 + plot2 + plot3
  dev.off()

  seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 50)
}


### RNA data with SCTransform
if (opt$batch == "null") {
  message(">>> Normalizing with SCTransform...")
  seu <- SCTransform(seu, verbose = TRUE)
} else {
  message(">>> Normalizing with SCTransform by batch...")
  seu.list <- SplitObject(seu, split.by = opt$batch)
  seu.list <- lapply(seu.list, SCTransform)
  lapply(seu.list, function(s) {
    # circumvent Seurat bug https://github.com/satijalab/seurat/pull/3658
    Misc(s[["SCT"]], "vst.out")$cells_step1 <- rownames(Misc(s[["SCT"]], "vst.out")$cell_attr)
  })
  message(">>> Integrating batches using SCTransform...")
  seu.features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 2000)
  seu.list <- PrepSCTIntegration(object.list = seu.list,
                                 anchor.features = seu.features,
                                 verbose = TRUE)
  seu.anchors <- FindIntegrationAnchors(object.list = seu.list,
                                        normalization.method = "SCT",
                                        anchor.features = seu.features,
                                        verbose = TRUE)
  seu <- IntegrateData(anchorset = seu.anchors,
                       normalization.method = "SCT",
                       verbose = FALSE)
  rm(seu.list)
  rm(seu.anchors)
  message(">>> Saving integrated data to ", debugdir)
  saveRDS(seu, file.path(debugdir, "integrated.rds"))
}


# If batch integration has occured, PCs, clusters, and UMAP will
# be caluclated based on integrated data.

message(">>> Running standard analyses on SCT data...")
message(">>> Calculating PCs...")
seu <- RunPCA(seu, verbose=TRUE)

### clustering based RNA
message(">>> finding clusters...")
seu <- FindNeighbors(seu, dims = 1:30, graph.name="rna.snn") %>%
  FindClusters(graph.name="rna.snn")
seu[["rnaClusterID"]] <- Idents(seu)

### RNA UMAP
message(">>> Calculating RNA UMAP")
seu <- RunUMAP(seu,
               reduction = 'pca',
               dims = 1:30,
               assay = 'RNA',
               reduction.name = 'rna.umap',
               reduction.key = 'rnaUMAP_')


# ----- Process ADT -----
adt <- "ADT" %in% Assays(seu)
if (adt) {
  message(">>> Processing ADT")
  DefaultAssay(seu) <- 'ADT'

  ### use all ADT features for dimensional reduction
  message(">>> Calculating CLR transformation...")
  VariableFeatures(seu) <- rownames(seu[["ADT"]])
  if (!opt$normalized) {
    seu <- NormalizeData(seu, normalization.method = 'CLR')
  }
  seu[["ADT"]]@data@x[is.na(seu[["ADT"]]@data@x)] <- 0
  seu[["ADT"]]@data@x[seu[["ADT"]]@data@x == -Inf] <- 0

  seu <- ScaleData(seu) %>% RunPCA(verbose=FALSE, reduction.name = 'adt.pca')

  ### clustering based ADT
  message(">>> finding clusters...")
  seu <- FindNeighbors(seu, dims = 1:20, graph.name="adt.snn") %>%
    FindClusters(graph.name="adt.snn")
  seu[["adtClusterID"]] <- Idents(seu)

  ### Identify multimodal neighbors
  message(">>> Calculating wnn graph...")
  seu <- FindMultiModalNeighbors(
    seu, reduction.list = list("pca", "adt.pca"),
    dims.list = list(1:30, 1:20), modality.weight.name = "RNA.weight"
  )

  ### clustering based on wnn
  message(">>> Finding wnn clusters...")
  seu <- FindClusters(seu, graph.name = "wsnn")
  seu[["wsnnClusterID"]] <- Idents(seu)

  ### UMAPs based on RNA, ADT, and WNN
  message(">>> Calculating UMAPs...")
  seu <- RunUMAP(seu, reduction = 'adt.pca', dims = 1:20, assay = 'ADT',
                 reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
  seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                 reduction.key = "wnnUMAP_")
}
message(">>> Saving processed data to ", debugdir)
saveRDS(seu, file.path(debugdir, "process.rds"))


# ----- Map to reference -----
message(">>> Mapping to reference dataset...")
DefaultAssay(seu) <- "SCT"

message(">>> Loading reference data from ", opt$referencedat)
reference <- LoadH5Seurat(file.path(opt$referencedat))

### find anchors between reference and query
### use precomputed supervised PCA (spca) transformation
message(">>> Finding anchors...")
anchors <- FindTransferAnchors(
  reference = reference,
  query = seu,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

message(">>> Adding mapping scores...")
seu <- AddMetaData(
  seu,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

message(">>> Mapping query dataset to reference...")
seu <- MapQuery(
  anchorset = anchors,
  query = seu,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# Add new UMAP
message(">>> Calculating new UMAP based on merged datasets...")
reference$id <- "reference"
seu$id <- "query"
refquery <- merge(reference, seu)
refquery[["spca"]] <- merge(reference[["spca"]], seu[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = "spca", dims = 1:50)
# Add new UMAP to misc slot in seurat object.
merged_umap <- refquery[["umap"]]
merged_umap@misc <- list(id = refquery$id)
seu@misc <- list(merged_umap = merged_umap)
rm(refquery)

# Save output
message(">>> Saving Seurat object...")
saveRDS(seu, file.path(outpath_seu))
message(">>> Saving h5Seurat object...")
SaveH5Seurat(seu, file.path(outpath_h5s), overwrite = TRUE)


# ----- Plots -----
### level1 cell type on rna, adt, wnn, and ref umaps
# pdf(file=file.path(plotdir, "celltypel1.pdf"),
#     paper="USr")
# p1 <- DimPlot(seu,
#               reduction = 'rna.umap',
#               group.by = 'predicted.celltype.l1',
#               label = TRUE,
#               repel = TRUE, label.size = 2.5) + NoLegend()
# if (adt) {
#   p2 <- DimPlot(seu,
#                 reduction = 'adt.umap',
#                 group.by = 'predicted.celltype.l1',
#                 label = TRUE,
#                 repel = TRUE, label.size = 2.5) + NoLegend()
#   p3 <- DimPlot(seu,
#                 reduction = 'wnn.umap',
#                 group.by = 'predicted.celltype.l1',
#                 label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
# } else {
#   p2 <- NULL
#   p3 <- NULL
# }
# p4 = DimPlot(seu,
#              reduction = "ref.umap",
#              group.by = "predicted.celltype.l1",
#              label = TRUE, label.size = 2.5, repel = TRUE) + NoLegend()
# p1+p2+p3+p4
# dev.off()
#
# ### l2 labels
# pdf(file=file.path(plotdir, "celltypel2.pdf"),
#     paper="USr")
# p1 <- DimPlot(seu,
#               reduction = 'rna.umap',
#               group.by = 'predicted.celltype.l2',
#               label = TRUE,
#               repel = TRUE, label.size = 2.5) + NoLegend()
# if (adt) {
#   p2 <- DimPlot(seu,
#                 reduction = 'adt.umap',
#                 group.by = 'predicted.celltype.l2',
#                 label = TRUE,
#                 repel = TRUE, label.size = 2.5) + NoLegend()
#   p3 <- DimPlot(seu,
#                 reduction = 'wnn.umap',
#                 group.by = 'predicted.celltype.l2',
#                 label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
# } else {
#   p2 <- NULL
#   p3 <- NULL
# }
#
# p4 = DimPlot(seu,
#              reduction = "ref.umap",
#              group.by = "predicted.celltype.l2",
#              label = TRUE, label.size = 2.5, repel = TRUE) + NoLegend()
# p1+p2+p3+p4
# dev.off()
#
# # Merged UMAP
# pdf(file=file.path(plotdir, "merged.pdf"),
#     paper="USr")
# ggplot(data.frame(merged_umap[[, 1:2]])) +
#   geom_point(aes(x = UMAP_1, y = UMAP_2, color = merged_umap@misc$id)) +
#   scale_color_discrete(name = "ID")
# dev.off()
