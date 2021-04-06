#!/bin/bash
ml fhR/4.0.3-foss-2020b
cd $wd/process
DATASET=$1

REFERENCEDAT=/fh/scratch/delete10/gottardo_r/hmiller/pbmc_multimodal.h5seurat
DATADIR=/fh/fast/gottardo_r/ytian_working/covid19_datasets
DEBUGDIR=/fh/scratch/delete10/gottardo_r/hmiller
PLOTDIR=/fh/fast/gottardo_r/ytian_working/covid19_datasets/plots_seurat

INDIR=${DATADIR}/seu
OUTDIR=${DATADIR}/seu_reprocessed
H5OUT=${DATADIR}/h5seurat_geneqc

Rscript process.R -D $DATASET -i $INDIR -o $OUTDIR -d $DEBUGDIR -r $REFERENCEDAT --plotdir $PLOTDIR --h5dir $H5OUT -b $batch
