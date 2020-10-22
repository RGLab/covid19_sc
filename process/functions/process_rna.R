process_rna <- function(seu) {
  ### basic QC
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 50)
  
  ### RNA data with SCTransform
  seu <- SCTransform(seu, verbose = FALSE)
  
  seu <- RunPCA(seu, verbose=FALSE)
  
  ### clustering based RNA
  seu <- FindNeighbors(seu, dims = 1:30, graph.name="rna.snn") %>%
    FindClusters(graph.name="rna.snn")
  seu[["rnaClusterID"]] <- Idents(seu)
  
  return(seu)
}