plot_results <- function(seu,
                         plotDir) {
  
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
}

