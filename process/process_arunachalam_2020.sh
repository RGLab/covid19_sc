#!/bin/bash
DATADIR=$1
REFERENCEDAT=$2

DATASET=arunachalam_2020
DEBUGDIR=${3}$DATASET
INDAT=${DATADIR}${DATASET}.step1.rds
OUTDAT=${DATADIR}${DATASET}_processed.rds
PLOTDIR=${4}${DATASET}
NCORES=${5}

Rscript process/process.R -i $INDAT -o $OUTDAT -d $DEBUGDIR -r $REFERENCEDAT --adt --plotdir $PLOTDIR