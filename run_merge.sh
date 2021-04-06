#!/bin/bash
#SBATCH --constraint=gizmok
#SBATCH --mail-user=hmiller@fredhutch.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/fh/fast/gottardo_r/ytian_working/covid19_datasets/covid19_sc/logs/merge-%J.log

ml fhR/4.0.3-foss-2020b
cd /fh/fast/gottardo_r/ytian_working/covid19_datasets/covid19_sc/merge
Rscript ./merge_datasets_new.R
Rscript ./merge_sce.R

ml fhPython/3.8.2-foss-2020a-Python-3.8.2
python merge_sct.py
