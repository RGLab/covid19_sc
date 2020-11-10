# covid19_sc
Processing workflow for covid19 single cell data

## renv
This project uses [renv](https://rstudio.github.io/renv/articles/renv.html) for package management. 
To activate, simply open an R session in the top-level directory of the repo. This will prompt the 
installation of the correct version of the `renv` package. Then call `renv::restore()` to install all packages to a project-specific library. 

If you wish to disable renv for this project, simply remove or comment out the following line from 
the project `.Rprofile`: 
```
source("renv/activate.R")
```

## Processing Methods
### 1. Process datasets. 
11 single cell datasets were processed using the Seurat package. Low-quality cells were filtered using basic QC thresholds, normalization was performed using SCtransform, a PCA was performed, clusters were identified using the leiden algorithm, and a UMAP was calculated for visualization. The results were saved to the final Seurat object. 

All datasets were mapped to a reference PBMC dataset, using [Seurat V4 methods](https://satijalab.org/seurat/v4.0/reference_mapping.html) and the PBMC reference supplied with Seurat V4. The output from this analysis is saved with the processed dataset, and includes the reference UMAP as `ref.umap`, predicted protein expression values as the assay `predicted_ADT`, and cell type annotations. 

Two datasets also included protein expression data. On these datasets, we also normalized protein expression data, ran a PCA on the result, and identified clusters from the PCA. Furthermore, a weighted nearest neighbor graph was calculated to integrate protein and gene expression data using Seurat V4, and clusters and a UMAP were calculated from the result. 

`process.R` runs basic QC and processing steps on a Seurat object. `process.sh` and `run.sh` are wrappers which will set off a slurm job to process a dataset. The 11 included datasets were processed using `run.sh` with the following options: 
  * `run.sh lee_2020 4`
  * `run.sh wilk_2020 8`
  * `run.sh chua_2020 16`
  * `run.sh meckiff_2020 16`
  * `run.sh arunachalam_2020 4`
  * `run.sh liao_2020 8`
  * `run.sh su_2020 16`
  * `run.sh yu_2020 16`
  * `run.sh wen_2020 4`
Two datasets were run with batch correction: 
  * `run.sh liao_2020 8 sample`
  * `run.sh zhu_2020 8 batch`
One dataset was too large to run the mapping step in one batch, so the mapping was done using the `map_to_reference.R` script: 
  ```
  run.sh su_2020 16
  map_to_reference_su_2020.sh
  ```
  
  
### 2. Standardize metadata. 
For each dataset, once `process.R` has completed and the processed version is saved, metadata is extracted from the resulting h5seurat file and written to a tsv file via `write_metadata.R`. Relevant standard fields were derived using the script `standardize_metadata.R`, which writes the output to a tsv file. Standard `disease_severity` scores and `days_since_symptom_onset` were assigned to each COVID-19 patient. For patients were either of these fields were missing, they were estimated as follows: 

#### `disease_severity`

The standardized severity score is based as closely as possible on this model, based 
on information provided in datasets and supplementary materials: 
1. `mild` = not hospitalized
2. `moderate` = hospitalized (no ICU)
3. `severe` = ICU

#### `days_since_symptom_onset`

Where provided, this is the number of days between symptom onset and the sampling
time. It was estimated as follows for patients where it was not provided: 
1. If days since hospitalization is known, 7 days are added to that number. This is based on the median 
difference between symptom onset and hospitalization for the 14 patients where both dates were reported.
2. For two datasets where neither days since hospitalization nor days since symptom
onset are known: 
  1. su_2020: `T1` (sample near the time of diagnosis) was assigned `7`, with the assumption that diagnosis
  and hospitalization would happen on a similar timeline. `T2` (approximately 
  one week later) was assigned `14`. 
  2. wen_2020: Early recovery stage (ERS) patients were assigned `15`. Late recovery stage (LRS) patients
  were assigned `29`. These numbers are rough estimates, based on the CDC guidelines for healthcare workers (https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html).
    1. The early-recovery stage (ERS) of COVID-19 patients are defined as the date of nucleic acid turning negative to blood sampling is less than seven days.
    2. The late-recovery stage (LRS) of COVID-19 patients are defined as the date of nucleic acid turning negative to blood sampling is more than fourteen days.


### 3. Merge datasets
Standard metadata, SCT expression values, and predicted protein scores were merged for all datasets. `merge_datasets.R` performs this task, and saves the result to an h5seurat file. 

Standard metadata was merged into each dataset using `merge_metadata.R`. 

