# Metadata standardization scripts

## Usage
These scripts must be run in the following order: 

1. `write_metadata.R`
  This script will pull metadata out of a processed dataset and write the results to a tsv file in `$DATA_DIR/metadata`. 
  Example: `Rscript write_metadata.R arunachalam_2020`

2. `standardize_metadata.R`
  This reads in the tsv files generated from `write_metadata.R` as well as any required supplementary files for each 
  dataset, then merges the results and writes it to `$DATA_DIR/metadata/all_meta.tsv`.
  
3. `merge_metadata.R`
  This will merge standard metadata back into the processed h5seurat file, for all datasets. Standard fields 
  get `_standard` appended to the field name. 
  
4. `add_column_attributes.R`
  This updates the metadata in the h5seurat files with column attributes indicating how they should be 
  displayed in the web portal, for all datasets.

---

`update_metadata.R` --  This is a helper script, which will merge updated standardized metadata into the _merged_ dataset. 
  
  
## Supplemental metadata
The files in the `supplemental` directory are pulled from supplemental files from the publications or provided by authors, and required for deriving standard fields. 
