#!/bin/bash

#' Cellranger Aggr
#'
#' Usage:
#' sbatch -n 1 -c 32 --mem=128G --time=120:00:00 process_10x_aggr.sh

## Base working directory
BASE="/fh/fast/gottardo_r/ytian_working/covid19_datasets/data/silvin_2020"

## load cellranger module
ml CellRanger/4.0.0

## aggr ---------------------------------------------------------
echo "Working on aggr..."

mkdir -p $BASE/cellranger/aggr
cd $BASE/cellranger/aggr

cellranger aggr \
	   --id=silvin_2020 \
	   --csv=$BASE/silvin_2020_aggr.csv 

echo "Run is complete."
