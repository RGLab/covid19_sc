process_adt <- function(seu) {
  DefaultAssay(seu) <- 'ADT'
  ### use all ADT features for dimensional reduction
  VariableFeatures(seu) <- rownames(seu[["ADT"]])
  seu <- NormalizeData(seu, normalization.method = 'CLR') %>% 
    ScaleData() %>% RunPCA(verbose=FALSE, reduction.name = 'adt.pca')
  
  ### clustering based ADT
  seu <- FindNeighbors(seu, dims = 1:20, graph.name="adt.snn") %>%
    FindClusters(graph.name="adt.snn")
  seu[["adtClusterID"]] <- Idents(seu)
  
  ### Identify multimodal neighbors
  seu <- FindMultiModalNeighbors(
    seu, reduction.list = list("pca", "adt.pca"), 
    dims.list = list(1:30, 1:20), modality.weight.name = "RNA.weight"
  )
  
  ### clustering based on wnn
  seu <- FindClusters(seu, graph.name = "wsnn")
  seu[["wsnnClusterID"]] <- Idents(seu)
  
  ### UMAPs based on RNA, ADT, and WNN
  seu <- RunUMAP(seu, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                 reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  seu <- RunUMAP(seu, reduction = 'adt.pca', dims = 1:20, assay = 'ADT', 
                 reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
  seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                 reduction.key = "wnnUMAP_")
  
  return(seu)
}