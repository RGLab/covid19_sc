# Scripts for merging datasets

Scripts in this directory were used for merging processed datasets. The final merged gene expression matrix is too large to be represented as a sparse matrix in R, so cannot be represented as a Seurat object. To merge all datasets:  

1. `merge_datasets.R` 
  This will create two Seurat objects with the standardized metadata and assays and save the result to h5seurat files. Some additional standardization occurs in this step:  
    * Gene names are standardized, using the most current mappings from the [HUGO Gene nomenclature committee](https://www.genenames.org/). This mapping table is provided in this directory as `hgncAlias2Symbol.tsv`  
    * Genes are filtered to only include those present in all datasets. For simplicity, this list of genes is provided in `genes.txt`.   

2. `merge_sce.R`  
  This merges the two Seurat objects created from `merge_datasets.R` utilizing the Bioconductor SingleCellExperiment object and HDF5Array, which allows matrices to be referenced from disk without reading into memory, and saves the result as a SingleCellExperiment object with the transcriptomic data saved in hdf5 format.   
  
3. `merge_sct.py`  
  This script merges the two Seurat objects created from `merge_datasets.R` utilizing the Python hdf5 API to coerce the merged object into the h5seurat format, which can be transformed into the web portal.  
  
  
