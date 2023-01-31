# sim-debiased-inference
Numerical studies for the paper "Debiased inference for a covariate-adjusted regression function"

# Reproduce the paper's numerical studies
This repository requires the R library [DebiasedDoseResponse](https://github.com/Kenta426/DebiasedDoseResponse). You can install this package by running the following prompt in R:
```
library(devtools)
devtools::install_github("Kenta426/DebiasedDoseResponse")
library(DebiasedDoseResponse)
```

Next, the following prompts run the exact experiments from the original paper:
```
# enter this prompt directly in Terminal
./sim_pointwise.sh  # Run simulation to generate Fig. 2, Fig. 3 on the arxiv version of the paper
./sim_effect.sh     # Run simulation to generate Fig. 4 on the arxiv version of the paper
./sim_uniform.sh    # Run simulation to generate Fig. 5 on the arxiv version of the paper
```
Note that the original results are based on 1000 simulations from seed `1000` to `2000` and ran in parallel 
on a remote machine and not recommended for a local machine. 

The simulation results will be stored in the following directories:
```
# the results from ./sim_pointwise.sh 
./result/pointwise/[n]/[seed]_[i].Rdata
# the results from ./sim_effect.sh 
./result/effect/[n]/[seed]_[i].Rdata
```
where [n] corresponds to the sample size, [seed] corresponds the base seed (i.e., one of {1000, 1100, ..., 1900}) and [i] corresponds the ith iteration of the seed (i.e., one of {1, 2, ..., 100}). For example, `./result/effect/2500/1200_35.Rdata` corresponds to the experiments from `./sim_effect.sh` when `n=2500` with seed `1234`.
Both `./sim_pointwise.sh` and `./sim_uniform.sh` store data in the same locations (i.e., the experiments only differ by the sample sizes). 

