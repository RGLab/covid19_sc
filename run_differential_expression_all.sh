# Starts an sbatch job for each dataset, which in turn starts sbatch jobs for
# each cell type to run differential expression analysis and identify the top
# 30 marker genes, writing them to the h5seurat file.

./run_differential_expression.sh arunachalam_2020 2
./run_differential_expression.sh bacher_2020 2
./run_differential_expression.sh chua_2020 4
./run_differential_expression.sh combes_2021 4
./run_differential_expression.sh grant_2021 4
./run_differential_expression.sh heming_2021 4
./run_differential_expression.sh kusnadi_2021 4
./run_differential_expression.sh lee_2020 4
./run_differential_expression.sh liao_2020 4
./run_differential_expression.sh meckiff_2020 4
./run_differential_expression.sh schulte-schrepping_2020 4
./run_differential_expression.sh silvin_2020 4
./run_differential_expression.sh stephenson_2021 8
./run_differential_expression.sh su_2020 8
./run_differential_expression.sh trump_2020 4
./run_differential_expression.sh wen_2020 4
./run_differential_expression.sh yao_2021 4
./run_differential_expression.sh yu_2020 4
./run_differential_expression.sh zhu_2020 4
./run_differential_expression.sh bost_2021 4
./run_differential_expression.sh combes_2021 4
