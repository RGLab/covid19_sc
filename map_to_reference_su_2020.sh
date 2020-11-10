#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G

INPUT=/fh/fast/gottardo_r/ytian_working/covid19_datasets/seu/su_2020_normalized.rds
OUTPUT=/fh/fast/gottardo_r/ytian_working/covid19_datasets/seu/su_2020_processed.rds
H5PATH=/fh/fast/gottardo_r/ytian_working/covid19_datasets/h5seurat/su_2020_processed.h5seurat

ml fhR/4.0.2-foss-2019b
cd /fh/fast/gottardo_r/ytian_working/covid19_datasets/covid19_sc
Rscript map_to_reference.R -i $INPUT -o $OUTPUT -s $H5PATH -b batch
