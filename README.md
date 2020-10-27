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
