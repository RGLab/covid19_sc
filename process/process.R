# ----- Parse options ------
library(optparse)
# TODO: add option for future plan
option_list = list(
  make_option(c("-i", "--inpath"),
              help = "path to load input seurat object",
              default = ""),
  make_option(c("-o", "--outpath"),
              help = "path to store final seurat object",
              default = ""),
  make_option(c("-d", "--debugdir"), 
              help = "path to store intermediate objects for debugging",
              default = ""),
  make_option(c("-r", "--referencedat"), 
              help = "path to reference data for mapping",
              default = ""),
  make_option(c("--plotdir"),
              help = "path to write plots",
              default = ""),
  make_option(c("-a", "--adt"),
              help = "process adt data",
              action = "store_true",
              default = TRUE)
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# check args
if (!file.exists(opt$inpath)) {
  stop("could not find ", opt$inpath)
} 
if (!dir.exists(opt$debugdir)) {
  message("creating debug dir: ", opt$debugdir)
  dir.create(opt$debugdir)
}
if (!dir.exists(opt$plotdir)) {
  dir.create(opt$plotdir)
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
  library(patchwork)
}) 

set.seed(1001)
ncores <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (ncores != "") {
  future::plan(strategy = "multicore", workers = as.integer(ncores), seed = 1001)
  # Set the maxSize for future to 2G per core
  options(future.globals.maxSize = 2 * as.integer(ncores) * 1024^3)
}

# ----- START -----
message(">>> Reading input: ", opt$inpath)
seu <- readRDS(opt$inpath)

### ----- Process RNA -----
message(">>> Processing RNA...")
message(">>> QC cells...")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 50)

### RNA data with SCTransform
message(">>> runing SCTransform...")
seu <- SCTransform(seu, verbose = FALSE)
message(">>> Calculating PCs...")
seu <- RunPCA(seu, verbose=FALSE)

### clustering based RNA
message(">>> finding clusters...")
seu <- FindNeighbors(seu, dims = 1:30, graph.name="rna.snn") %>%
  FindClusters(graph.name="rna.snn")
seu[["rnaClusterID"]] <- Idents(seu)


# ----- Process ADT -----
if (opt$adt) {
  message(">>> Processing ADT")
  DefaultAssay(seu) <- 'ADT'
  
  ### use all ADT features for dimensional reduction
  message(">>> Calculating CLR transformation...")
  VariableFeatures(seu) <- rownames(seu[["ADT"]])
  seu <- NormalizeData(seu, normalization.method = 'CLR') %>% 
    ScaleData() %>% RunPCA(verbose=FALSE, reduction.name = 'adt.pca')
  

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
  seu <- RunUMAP(seu, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                 reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  seu <- RunUMAP(seu, reduction = 'adt.pca', dims = 1:20, assay = 'ADT', 
                 reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
  seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                 reduction.key = "wnnUMAP_")
}
message(">>> Saving processed data to ", opt$debugdir)
saveRDS(seu, file.path(opt$debugdir, "process.rds"))


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

message(">>> Mapping query dataset to reference...")
seu <- MapQuery(
  anchorset = anchors,
  query = seu,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

# Add new UMAP
message(">>> Calculating new UMAP based on merged datasets...")
reference$id <- 'reference'
seu$id <- 'query'
refquery <- merge(reference, seu)
refquery[["spca"]] <- merge(reference[["spca"]], seu[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
# Add new UMAP to misc slot in seurat object. 
merged_umap <- refquery[["umap"]]
merged_umap@misc <- list(id = refquery$id)
seu@misc <- list(merged_umap = merged_umap)
rm(refquery)

saveRDS(seu, file.path(opt$outpath))


# ----- Plots -----
plotDir <- opt$plotdir
### level1 cell type on rna, adt, wnn, and ref umaps
pdf(file=file.path(plotDir, "celltypel1.pdf"),
    paper="USr")
p1 <- DimPlot(seu,
              reduction = 'rna.umap',
              group.by = 'predicted.celltype.l1',
              label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(seu,
              reduction = 'adt.umap',
              group.by = 'predicted.celltype.l1',
              label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 <- DimPlot(seu,
              reduction = 'wnn.umap',
              group.by = 'predicted.celltype.l1',
              label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p4 = DimPlot(seu,
             reduction = "ref.umap",
             group.by = "predicted.celltype.l1",
             label = TRUE, label.size = 2.5, repel = TRUE) + NoLegend()
p1+p2+p3+p4
dev.off()

### l2 labels
pdf(file=file.path(plotDir, "celltypel2.pdf"),
    paper="USr")
p1 <- DimPlot(seu,
              reduction = 'rna.umap',
              group.by = 'predicted.celltype.l2',
              label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(seu,
              reduction = 'adt.umap',
              group.by = 'predicted.celltype.l2',
              label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 <- DimPlot(seu,
              reduction = 'wnn.umap',
              group.by = 'predicted.celltype.l2',
              label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p4 = DimPlot(seu,
             reduction = "ref.umap",
             group.by = "predicted.celltype.l2",
             label = TRUE, label.size = 2.5, repel = TRUE) + NoLegend()
p1+p2+p3+p4
dev.off()


