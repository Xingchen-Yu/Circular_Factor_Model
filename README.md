# This repository is superseded by https://github.com/Xingchen-Yu/SLFM1D. 

# Circular Factor Model
[![DOI](https://zenodo.org/badge/271097395.svg)](https://zenodo.org/badge/latestdoi/271097395)
## Data

Roll call data for the U.S. House of Representatives (100th  to 116th) is obtained from voteview.com.  Each data is binary and contains some missing values. 

## Code
Circular_Factor_Model.R takes 4 arguments which are iterations, burn-in period of mcmc, number of cores used, the number of the U.S
House of representatives data. It is recommended to run on a server with many cores.

For exammple, to run 100000 iterations with 80000 burn-in using 12 cores for the 112th U.S House of representative data,
one could excute the following command.

Rscript ./main_script/Circular_Factor_Model.R 100000 80000 12 112 &

The “main_script” folder contains the following scripts,

Circular_Factor_Model.R implements the proposed model in the paper. (Parallel computing required)

Euclidean_1d. R implements the 1d Eulicdean latent factor model in the paper. (Single thread)

Euclidean_2d. R implements the 2d Eulicdean latent factor model in the paper. (Single thread)

tables_figures.R reproduces the major results in the paper once the corresponding runs are finished

The “source” folder contains the following scripts,

circular_factor_model_functions. R includes cpp functions for the proposed model.

read_kh2.R and ymat_spit.R contains data preprocessing and two miscellaneous functions

## Reproducibility workflow
Under the root directory, first run each the following commands,

### Run circular factor model for the 112 House data, with 100000 iterations, 80000 burn-in and 12 cores (2-3 days)

Rscript ./main_script/Circular_Factor_Model.R 100000 80000 12 112 &

### Run circular factor model for the 116 House data, with 100000 iterations, 80000 burn-in and 12 cores (1-2 days)

Rscript ./main_script/Circular_Factor_Model.R 100000 80000 12 116 &

### Run 1d Euclidean latent factor model for the 112 House data with 30000 iterations and 10000 burnin

Rscript ./main_script/Euclidean_1d.R 30000 10000 112 &

### Run 1d Euclidean latent factor model for the 116 House data with 30000 iterations and 10000 burnin

Rscript ./main_script/Euclidean_1d.R 30000 10000 116 &

### Run 2d Euclidean latent factor model for the 116 House data with 40000 iterations and 20000 burnin

Rscript ./main_script/Euclidean_2d.R 40000 20000 116 &

There will be a total of 7 Rdata files stored in the same directory. Tables 1, 2, 3 and Figures 3, 4, 5 can then be generated by running the tables_figures.R script.

## System and pacakge requirements
### Versions of softwares
R version 3.6.0 (2019-04-26)

Platform: x86_64-redhat-linux-gnu (64-bit)

Running under: CentOS Linux 7 (Core)

### Libraries and dependencies
wnominate_1.2.5, pscl_1.5.5, truncnorm_1.0-8, mvnfast_0.2.5
RcppArmadillo_0.9.300.2.0, rlecuyer_0.3-5, snowfall_1.84-6.1, 
snow_0.4-3,   Rcpp_1.0.1 (or later)  






