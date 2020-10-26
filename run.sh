#!/bin/bash
# Usage: ./run.sh dataset_name ncores
# if ncores is not specified, defaults to 1. 

sbatch --constraint=gizmok \
  -J ${1} \
  -o ${0%/*}/logs/${1}_$(date +"%Y-%m-%d-%H:%M:%S").log \
  --mail-user=hmiller@fredhutch.org \
  --mail-type=ALL \
  --cpus-per-task=${2-1} \
  --time=08:00:00 \
  --export=wd=${0%/*} \
  process.sh ${1} 
  
