#!/bin/bash
ml fhR/4.0.2-foss-2019b
cd $wd
DATASET=$1

REFERENCEDAT=/fh/fast/gottardo_r/ytian_working/covid19_datasets/data/seurat/pbmc_multimodal.h5seurat
DATADIR=/fh/fast/gottardo_r/ytian_working/covid19_datasets
DEBUGDIR=/fh/scratch/delete10/gottardo_r/hmiller
PLOTDIR=/fh/fast/gottardo_r/ytian_working/covid19_datasets/plots_seurat/$DATASET

INDAT=${DATADIR}/seu/${DATASET}.step1.rds
OUTDAT=${DATADIR}/seu/${DATASET}_processed.rds
H5OUT=${DATADIR}/h5seurat/${DATASET}_processed.h5seurat

Rscript process.R -i $INDAT -o $OUTDAT -d $DEBUGDIR -r $REFERENCEDAT --plotdir $PLOTDIR --hpath $H5OUT
