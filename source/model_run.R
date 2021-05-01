#### checking and installing required packages###
required_package = c('Rcpp','snowfall','tidyverse','rlecuyer','RcppArmadillo','matrixStats')
check_package = sum(unlist(lapply(required_package, require, character.only = TRUE)))== length(required_package)
if(check_package ==F){
  install.packages(required_package,repos = "http://cran.us.r-project.org")
  lapply(required_package, require, character.only = TRUE)
}


source(file="./source/tidyverse_load_data.R")
source(file='./main_script/Circular_Factor_Model_tidyverse.R')
################################################################
hn = 116
house = T
h_s = ifelse(house==T,'H','S')
out = get_rollcall_data(house,h_s,hn,threshold = 0.4)

################################################################
# ymat = as.matrix(out[[1]] %>% select(-bioname,-name_district,-icpsr))
# pol_info = out[[1]] %>% select(bioname,name_district,icpsr)
# dup_name = out[[2]]
# filtered_legislator = out[[3]]

master = SLFM(out, n_pos=20000,burnin=80000,thin = 1, hyperparams=list(a = 1, b = 1/10, ccc_a = 1, ccc_b=25, kappa_a = 1, omega_sd=0.1, kappa_sd=0.5,
                                                                    i_epi_lower = 0.005, i_epi_upper = 0.04, j_epi_lower = 0.01 ,j_epi_upper = 0.105,
                                                                    i_leap = 10, j_leap = 10,skip = 50, jitter = T, WAIC_group = T),
                                                                    initial_values=NULL,core=10,cluster_seed=8888)

str(master)
