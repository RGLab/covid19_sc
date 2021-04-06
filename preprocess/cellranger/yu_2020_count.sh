#!/bin/bash

#' Cellranger Count Processing
#'
#' Usage:
#' sbatch -n 1 -c 32 -p campus-new --mem=128G --time=120:00:00 process_10x_count.sh
#'
#' Triggers processing via cellranger of 10X data, from custom reference
#' creation to counting to aggregation. Custom for the project to figure
#' out parameters for automation of a pipeline.

## Base working directory
BASE="/fh/fast/gottardo_r/ytian_working/covid19_datasets/data/yu_2020"
FASTQ_DIR="${BASE}/fastq_clean"
REFERENCE_GENOME="/fh/fast/gottardo_r/ytian_working/run_cellranger_count/refdata-gex-GRCh38-2020-A"

## load cellranger module
ml CellRanger/4.0.0

## Count Quantification ---------------------------------------------------------
echo "Working on counts quantification..."

## Count in each sample
mkdir -p $BASE/cellranger/count
cd $BASE/cellranger/count

for SAMPLE_DIR in $(ls -d ${FASTQ_DIR}/*/); do
	echo "${SAMPLE_DIR}"
	cellranger count \
		--id=$(basename ${SAMPLE_DIR}) \
          	--transcriptome=${REFERENCE_GENOME} \
          	--fastqs=${SAMPLE_DIR} \
          	--localmem=128 \
	  	--localcores=32
done

echo "Run is complete."
