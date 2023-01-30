#!/bin/bash
#
echo "Start experiment"
for seed in 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900
do
    Rscript pointwise-inference.R 250 $seed
    Rscript pointwise-inference.R 500 $seed
    Rscript pointwise-inference.R 750 $seed
    Rscript pointwise-inference.R 1000 $seed
    Rscript pointwise-inference.R 1250 $seed
    Rscript pointwise-inference.R 1500 $seed
    Rscript pointwise-inference.R 2000 $seed
done
