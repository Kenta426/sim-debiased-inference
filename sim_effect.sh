#!/bin/bash
#
echo "Start experiment"
for seed in 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900
do
    Rscript effect-inference.R 500 $seed
    Rscript effect-inference.R 1000 $seed
    Rscript effect-inference.R 2500 $seed
done
