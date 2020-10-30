#!/bin/bash
# Usage: ./run.sh dataset_name ncores batch
# batch is optional

sbatch --constraint=gizmok \
  -J ${1} \
  -o ${0%/*}/logs/${1}_$(date +"%Y-%m-%d-%H:%M:%S").log \
  --mail-user=$(whoami)@fredhutch.org \
  --mail-type=ALL \
  --cpus-per-task=$2 \
  --time=48:00:00 \
  --export=wd=${0%/*},batch=${3-null}\
  process.sh ${1}

