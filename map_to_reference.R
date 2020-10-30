# sbatch --constraint=gizmok -J su_2020 -o logs/su_2020_map.log --mail-user=$(whoami)@fredhutch.org --mail-type=ALL --time=08:00:00 map_to_reference_su_2020.sh

# ----- Parse options ------
message(getwd())
library(optparse)
option_list = list(
  make_option(c("-i", "--inpath"),
              help = "path to load input seurat object",
              default = "/fh/scratch/delete10/gottardo_r/hmiller/liao_2020/process.rds"),
  make_option(c("-o", "--outpath"),
              help = "path to store final seurat object",
              default = "/fh/fast/gottardo_r/ytian_working/covid19_datasets/seu/liao_processed.rds"),
  make_option(c("-r", "--referencedat"),
              help = "path to reference data for mapping",
              default = "/fh/fast/gottardo_r/ytian_working/covid19_datasets/data/seurat/pbmc_multimodal.h5seurat"),
  make_option(c("-s", "--hpath"),
              help = "path to save h5Seurat object",
              default = "/fh/fast/gottardo_r/ytian_working/covid19_datasets/h5seurat/liao_2020_processed.h5seurat"),
  make_option(c("-b", "--batch"),
              help = "(optional) batch field. If specified, batches will be aligned using SCTransform methods.",
              default = "null")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

message(paste0(capture.output(opt), collaps = "\n"))

# ----- Map to reference -----
library(Seurat)
library(SeuratDisk)

message(">>> Mapping to reference dataset...")
seu <- readRDS(opt$inpath)
DefaultAssay(seu) <- "SCT"

message(">>> Loading reference data from ", opt$referencedat)
reference <- LoadH5Seurat(file.path(opt$referencedat))

if (opt$batch != "null") {
  seuList <- SplitObject(seu, opt$batch)
  anchors <- list()
  for (i in 1:length(seuList)) {
    message(">>> Finding anchors for batch ", i)
    anchors[[i]] <- FindTransferAnchors(
      reference = reference,
      query = seuList[[i]],
      normalization.method = "SCT",
      reference.reduction = "spca",
      dims = 1:50
    )
  }
  for (i in 1:length(seuList)) {
    message(">>> Mapping batch ", i)
    seuList[[i]] <- MapQuery(
      anchorset = anchors[[i]],
      query = seuList[[i]],
      reference = reference,
      refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2",
        predicted_ADT = "ADT"
      ),
      reference.reduction = "spca",
      reduction.model = "wnn.umap"
    )
  }
    message(">>> Saving seuList to /fh/scratch/delete10/gottardo_r/hmiller/su_2020/seuList_mapped.rds")
    saveRDS(seuList, "/fh/scratch/delete10/gottardo_r/hmiller/su_2020/seuList_mapped.rds")
  seu <- merge(seuList[[1]], seuList[2:length(seuList)], merge.dr = "ref.umap")
} else {

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
}

# Save output
message(">>> Saving h5Seurat object...")
SaveH5Seurat(seu, file.path(opt$hpath), overwrite = TRUE)

