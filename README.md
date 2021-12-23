# covid19_sc
Processing workflow for covid19 single cell data 
* Publication: [Single-cell immunology of SARS-CoV-2 infection](https://www.nature.com/articles/s41587-021-01131-y)
* Explore and download data: [Atlas](https://atlas.fredhutch.org/fredhutch/covid/)

## Setup
### renv
This project uses [renv](https://rstudio.github.io/renv/articles/renv.html) for package management. 
To activate, simply open an R session in the top-level directory of the repo. This will prompt the 
installation of the correct version of the `renv` package. Then call `renv::restore()` to install all packages to a project-specific library. 

If you wish to disable renv for this project, simply remove or comment out the following line from 
the project `.Rprofile`: 
```
source("renv/activate.R")
```
Scripts are designed to be run from the project directory. This will enable `renv` as well as load environment variables from `.Renviron`. 

### Paths and directory structure
Paths to data directories are set in the project-local `.Renviron`. All R scripts set paths relative to paths supplied in `.Renviron`. If you wish to run any scripts on your machine, simply edit the `.Renviron` file to point to the appropriate directories. There are also several bash scripts which wrap R scripts. Set paths manually in these scripts if you wish to run them. `.Renviron` variables are used as follows: 

* `DATA_DIR` -- All objects are saved and loaded in locations relative to this directory. It follows this format:  
```
DATA_DIR
├── data/
│  └── __raw data for all datasets__
├── h5seurat_reprocessed/
│  └── __output h5seurat files__
├── metadata/
│  └── __metadata is written here for intermediate standardization__
├── plots
│  └── __process.R writes some diagnostic plots here__
└── seu_reprocessed/
   └── __Preprocessed and processed Seurat objects are written here as .rds files__
```

* `REFERENCE_DATASET` -- Path to a local version of the reference dataset from [Hao et al](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1) available for download [here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat). 

* `DEBUG_DIR` -- Path to a directory for saving intermediate versions of the data for debugging. Utilized by `process.R`. 

## Processing 

The outline for processing all datasets are as follows: 

1. Pull datasets (sources noted below) and convert into Seurat objects using scripts in the `preprocess/` directory. 
1. `run_all.sh`: Process counts data for all datasets. Perform processing steps on assay data, including normalization, dimension reduction, cluster identification, cell type annotation, and reference mapping. More details below. 
1. Add cell type annotations to whole blood datasets using singleR. This is used when creating the combined dataset, as non-PBMC cells are removed. 
1. `run_update_metadata.sh`: Standardize metadata and update processed objects. More details below and in the `metadata/` directory. 
1. `run_merge.sh`: Merge PBMC and whole blood datasets into a combined dataset. More details below and in the `merge/` directory. 
1. `run_differential_expression_all.sh`: Identify differentially expressed marker genes for all cell types. Results are saved as csv files and merged into the `misc` slot in the h5seurat objects. 


### 1. Preprocess
First we downloaded data from public sources. Links for each source dataset are included below. Source data was prepared for processing by loading into a Seurat object and saving the result as an rds file. Scripts for preprocessing each dataset are provided in the `preprocess/` directory. 

### 2. Process datasets. 
21 single cell datasets were processed using the Seurat package. Low-quality cells were filtered using basic QC thresholds, normalization was performed using SCtransform, a PCA was performed, clusters were identified using the leiden algorithm, and a UMAP was calculated for visualization. The results were saved to the final Seurat object. 

All datasets were mapped to a reference PBMC dataset, using [Seurat V4 methods](https://satijalab.org/seurat/v4.0/reference_mapping.html) and the PBMC reference supplied with Seurat V4. The output from this analysis is saved with the processed dataset, and includes the reference UMAP as `ref.umap`, predicted protein expression values as the assay `predicted_ADT`, and cell type annotations. 

Three datasets also included protein expression data. On these datasets, we also normalized protein expression data, ran a PCA on the result, and identified clusters from the PCA. Furthermore, a weighted nearest neighbor graph was calculated to integrate protein and gene expression data using Seurat V4, and clusters and a UMAP were calculated from the result. 

`process.R` runs basic QC and processing steps on a Seurat object. `process.sh` and `run.sh` are wrappers which will set off a slurm job to process a dataset. `run_all.sh` will set off slurm jobs for all datasets, and can be used as a reference for which parameters were used for different datasets.


### 3. Standardize metadata. 
For each dataset, once `process.R` has completed and the processed version is saved, metadata is extracted and standardized, with standard fields merged back into the original h5seurat file. For more details about how the metadata scripts work, see the [README.md](/metadata/README.md) file in the metadata directory. The bash `run_update_metadata.sh` will run all the metadata standardization scripts in order, and can be used as a reference for how metadata standardization was done.

Standard `disease_severity` scores and `days_since_symptom_onset` were assigned to each COVID-19 patient. For patients were either of these fields were missing, they were estimated as follows: 

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
onset are known for one or more COVID-19 patients: 
  1. su_2020: `T1` (sample near the time of diagnosis) was assigned `7`, with the assumption that diagnosis
  and hospitalization would happen on a similar timeline. `T2` (approximately 
  one week later) was assigned `14`. 
  1. wen_2020: Early recovery stage (ERS) patients were assigned `15`. Late recovery stage (LRS) patients
  were assigned `29`. These numbers are rough estimates, based on the CDC guidelines for healthcare workers (https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html).  
    1. The early-recovery stage (ERS) of COVID-19 patients are defined as the date of nucleic acid turning negative to blood sampling is less than seven days.  
    2. The late-recovery stage (LRS) of COVID-19 patients are defined as the date of nucleic acid turning negative to blood sampling is more than fourteen days.  
  1. yao_2021: Samples from the same severity group were pooled, and due to an issue with hashing, while severity group each for each sample is known, they cannot be mapped to specific patients. Thus, we used the median days since symptom onset and hospitalization for each severity group.   
  1. bost_2021:  Authors reported the time from symptom onset to admission, but not the time until sampling. We assumed that samples were taken on the same day as admission, so used the median of the number of days between symptom onset and admission for each severity group, reported in Table 2.  
  1. stephenson_2021 and combes_2021: Both datasets were missing information about days since symptom onset for several patients. For those missing records, we used the median from the corresponding severity group where those numbers were not missing.  
    
    
### 4. Merge datasets
As only cell types in the Seurat v4 multimodal CITE-seq PBMC reference set can be correctly mapped, for visualization purposes we next took only the PBMC and whole blood datasets (16 in total, the latter minus neutrophils and basophils) and created a “merged dataset” consisting of over 2.5 million cells. This includese standard metadata, UMAP, SCT expression values, and predicted protein scores. 

The full, merged counts matrix is too large to be stored as a sparse matrix in R due to technical limitations. Thus, there are several steps required to create the hdf5 file required for the website. More details are available in the [README.md](/merge/README.md) file in the `merge` directory. The file `merge.sh` performs all the steps. 

## Data availability
Links for data sources as provided by the original publication are available below. Processed datasets ready for analysis are available at https://atlas.fredhutch.org/fredhutch/covid/ as Seurat objects or [h5seurat](https://mojaveazure.github.io/seurat-disk/articles/h5Seurat-load.html) files. 

### Data sources
Data were downloaded from the following locations: 

|dataset name            |data source                                                                                                                                                |
|:-----------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------|
|wilk_2020               |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728                                                                                               |
|yao_2021                |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154567                                                                                               |
|grant_2021              |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155249                                                                                               |
|chua_2020               |https://figshare.com/articles/COVID-19_severity_correlates_with_airway_epithelium-immune_cell_interactions_identified_by_single-cell_analysis/12436517     |
|bost_2021               |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157344                                                                                               |
|yu_2020                 |https://bigd.big.ac.cn/gsa/browse/CRA002572                                                                                                                |
|silvin_2020             |https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9221/                                                                                                |
|combes_2021             |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163668                                                                                               |
|trump_2020              |https://figshare.com/articles/dataset/Hypertension_delays_viral_clearance_and_exacerbates_airway_hyperinflammation_in_patients_with_COVID-19/13200278      |
|meckiff_2020            |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152522                                                                                               |
|wen_2020                |https://bigd.big.ac.cn/bioproject/browse/PRJCA002413                                                                                                       |
|lee_2020                |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689                                                                                               |
|bacher_2020             |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162086                                                                                               |
|su_2020                 |https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9357/                                                                                                |
|heming_2021             |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163005                                                                                               |
|schulte-schrepping_2020 |https://ega-archive.org/studies/EGAS00001004571                                                                                                            |
|kusnadi_2021            |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153931                                                                                               |
|liao_2020               |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926                                                                                               |
|zhu_2020                |https://db.cngb.org/search/project/CNP0001102/                                                                                                             |
|arunachalam_2020        |https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155673                                                                                               |
|stephenson_2021         |https://www.covid19cellatlas.org/index.patient.html                                                                                                        |



We acknowledge the Scientific Computing Infrastructure at Fred Hutch funded by ORIP grant S10OD028685.
