#!/bin/bash
#SBATCH --constraint=gizmok
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=/fh/fast/gottardo_r/ytian_working/covid19_datasets/covid19_sc/logs/metadata-%J.log

ml fhR/4.0.3-foss-2020b
cd /fh/fast/gottardo_r/ytian_working/covid19_datasets/covid19_sc/metadata

# Write metadata to metadata dir
Rscript write_metadata.R lee_2020
Rscript write_metadata.R wilk_2020
Rscript write_metadata.R chua_2020
Rscript write_metadata.R meckiff_2020
Rscript write_metadata.R arunachalam_2020
Rscript write_metadata.R silvin_2020
Rscript write_metadata.R su_2020
Rscript write_metadata.R yu_2020
Rscript write_metadata.R wen_2020
Rscript write_metadata.R liao_2020
Rscript write_metadata.R zhu_2020
Rscript write_metadata.R schulte-schrepping_2020
Rscript write_metadata.R bacher_2020
Rscript write_metadata.R kusnadi_2021
Rscript write_metadata.R heming_2021
Rscript write_metadata.R trump_2020
Rscript write_metadata.R yao_2021
Rscript write_metadata.R combes_2021
Rscript write_metadata.R grant_2021
Rscript write_metadata.R stephenson_2021

# Standardize metadata. This will write all_meta.tsv
Rscript standardize_metadata.R

# merge metadata back into processed h5seurat file, with _standard suffix
Rscript merge_metadata.R

# Add column attributes to h5seurat file indicating display type for portal
Rscript add_column_attributes.R
