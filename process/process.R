# ----- Parse options ------
library(optparse)
option_list = list(
  make_option(c("-i", "--inpath"),
              help = "path to load input seurat object",
              default = "/fh/fast/gottardo_r/ytian_working/covid19_datasets/seu/arunachalam_2020.step1.rds"),
  make_option(c("-o", "--outpath"),
              help = "path to store final seurat object"),
  make_option(c("-d", "--debugdir"), 
              help = "path to store intermediate objects for debugging",
              default = "/fh/scratch/delete10/gottardo_r/hmiller/arunachalam_2020"),
  make_option(c("-r", "--referencedat"), 
              help = "path to reference data for mapping",
              default = "/fh/fast/gottardo_r/ytian_working/covid19_datasets/data/seurat/pbmc_multimodal.h5seurat"),
  make_option(c("-a", "--adt"),
              help = "process adt data",
              action = "store_true",
              type = "logical", 
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
if (!file.exists(opt$referencedat)) {
  stop("could not find ", opt$referencedat)
}

# ----- Start -----

# Load functions and libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(magrittr)
  library(SeuratDisk)
  library(ggplot2)
  library(patchwork)
}) 
fundir <- "process/functions"
funScripts <- list.files(fundir)
for (file in funScripts) {
  source(file.path(fundir, file))
}

# Load raw seu object
seu <- readRDS(opt$inpath)

seu <- process_rna(seu)
if (opt$adt) {
  seu <- process_adt(seu)
}
seu <- map_to_reference(seu)
seu <- plot_results(seu)
