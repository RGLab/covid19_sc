# covid19_sc
Processing workflow for covid19 single cell data

## renv
This project uses [renv](https://rstudio.github.io/renv/articles/renv.html) for package management. 
To activate, simply open an R session in the top-level directory of the repo. This will prompt the 
installation of the correct version of the `renv` package. Then call `renv::restore()` to install 
all packages to a project-specific library. 

If you wish to disable renv for this project, simply remove or comment out the following line from 
the project `.Rprofile`: 
```
source("renv/activate.R")
```

# Standardized, merged dataset
On merging all datasets, we made several decisions to standardize metadata, as 
outlined below. 

### `disease_severity`

The standardized severity score is based as closely as possible on this model, based 
on information provided in datasets and supplementary materials: 
1. `mild` = not hospitalized
2. `moderate` = hospitalized (no ICU)
3. `severe` = ICU

### `days_since_symptom_onset`

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


## Processing datasets
`process.R` will process a dataset. For details: 
```
Rscript process.R --help
```

To start a slurm job: 
```
./run.sh <dataset_name>
```

## Processing notes
Batches: 
liao_2020 8 sample  
zhu_2020 8 batch
