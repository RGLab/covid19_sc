map_to_reference <- function(seu,
                             path_to_reference,
                             merge_and_reproject = FALSE) {
  DefaultAssay(seu) <- "SCT"
  
  reference <- LoadH5Seurat(file.path(path_to_reference))
  
  ### find anchors between reference and query
  ### use precomputed supervised PCA (spca) transformation
  anchors <- FindTransferAnchors(
    reference = reference,
    query = seu,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  
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
  
  if (merge_and_reproject) {
    reference$id <- 'reference'
    seu$id <- 'query'
    refquery <- merge(reference, seu)
    refquery[["spca"]] <- merge(reference[["spca"]], seu[["ref.spca"]])
    refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
    DimPlot(refquery, group.by = 'id')
    return(refquery)
  }
  
  return(seu)
}