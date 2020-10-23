#!/bin/bash
wd=${0%/*}
ml fhR/4.0.2-foss-2019b
cd $wd

DATASET=$1
REFERENCEDAT=/fh/fast/gottardo_r/ytian_working/covid19_datasets/data/seurat/pbmc_multimodal.h5seurat
DATADIR=/fh/fast/gottardo_r/ytian_working/covid19_datasets/seu/
DEBUGDIR=/fh/scratch/delete10/gottardo_r/hmiller/
PLOTDIR=/fh/fast/gottardo_r/ytian_working/covid19_datasets/plots_seurat/

$wd/process/process_${DATASET}.sh $DATADIR $REFERENCEDAT $DEBUGDIR $PLOTDIR 