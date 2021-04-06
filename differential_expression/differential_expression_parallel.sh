#!/bin/bash

ml fhR/4.0.3-foss-2020bml fhR/4.0.3-foss-2020b
cd $wd/differential_expression
Rscript differential_expression/differential_expression_parallel.R -D $1 ${2-4}
