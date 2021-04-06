#!/bin/bash
# Will start a job to run differential expression
# on a study, which will in turn start separate jobs
# for each cell type and run them in parallel.
# Specify dataset and number of cores per job (spawned by the main job)
# Usage: run_differential_expression.sh <dataset name> <cores per job>
# example: run_differential_expression.sh arunachalam_2020 2

# TODO: Check wd

sbatch --constraint=gizmok \
  -J DE_${1} \
  --export=wd=${0%/*} \
  -o ${0%/*}/logs/differential_expression_analysis/DE_${1}_$(date +"%Y-%m-%d-%H:%M:%S").log \
  --mail-user=$(whoami)@fredhutch.org \
  --mail-type=ALL \
  --cpus-per-task=1 \
  --time=48:00:00 \
  ./differential_expression/differential_expression_parallel.sh ${1} ${2-4}
