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
BASE="/fh/fast/gottardo_r/ytian_working/covid19_datasets/data/silvin_2020"
FASTQ_DIR="${BASE}/E-MTAB-9221"
REFERENCE_GENOME="/fh/fast/gottardo_r/ytian_working/run_cellranger_count/refdata-gex-GRCh38-2020-A"

SAMPLE_LABELS=(UPN2 UPN3 UPN5 UPN4 UPN4_j10 UPN1 UPN1_j10 UPN6 UPN6_j10 UPN7_j17)

## load cellranger module
ml CellRanger/4.0.0

## Count Quantification ---------------------------------------------------------
echo "Working on counts quantification..."

## Count in each sample
mkdir -p $BASE/cellranger/count
cd $BASE/cellranger/count

for I in $(seq 0 $((${#SAMPLE_LABELS[@]} - 1))); do
	echo "${SAMPLE_LABELS[$I]}"
	cellranger count \
		--id=${SAMPLE_LABELS[$I]} \
          	--transcriptome=${REFERENCE_GENOME} \
          	--fastqs=${FASTQ_DIR} \
		--sample=${SAMPLE_LABELS[$I]} \
          	--localmem=128 \
	  	--localcores=32
done

echo "Run is complete."
