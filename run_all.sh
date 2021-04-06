#!/bin/bash
# Usage: ./run_all.sh
# Sets off batch jobs for ALL datasets


./run.sh lee_2020 4
./run.sh wilk_2020 8
./run.sh chua_2020 16
./run.sh meckiff_2020 16
./run.sh arunachalam_2020 4
./run.sh silvin_2020 4
./run.sh su_2020 16
./run.sh yu_2020 16
./run.sh wen_2020 4

./run.sh liao_2020 8 sample
./run.sh zhu_2020 8 batch

./run.sh schulte-schrepping_2020 8 cohort
./run.sh bacher_2020 4
./run.sh kusnadi_2021 4
./run.sh heming_2021 8
./run.sh trump_2020 8
./run.sh yao_2021 4
./run.sh combes_2021 4

# These two were already normalized
./run_norm.sh grant_2021 8
./run_norm.sh stephenson_2021 16
