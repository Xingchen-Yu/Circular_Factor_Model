#### checking and installing required packages###
######################
required_package = c('Rcpp','snowfall','tidyverse','rlecuyer','RcppArmadillo','matrixStats')
check_package = sum(unlist(lapply(required_package, require, character.only = TRUE)))== length(required_package)
if(check_package ==F){
  install.packages(required_package,repos = "http://cran.us.r-project.org")
  lapply(required_package, require, character.only = TRUE)
}


source(file="./source/tidyverse_load_data.R")
source(file='./main_script/Circular_Factor_Model_tidyverse.R')
################################################################
hn = 117
house = T
h_s = ifelse(house==T,'H','S')
out = get_rollcall_data(house,h_s,hn,threshold = 0.4)

################################################################
# ymat = as.matrix(out[[1]] %>% select(-bioname,-name_district,-icpsr))
# pol_info = out[[1]] %>% select(bioname,name_district,icpsr)
# dup_name = out[[2]]
# filtered_legislator = out[[3]]
set.seed(2021)
out = matrix(sample(c(0,1),1000, T, prob=c(0.7,0.3)),100,10)
master = SLFM(out, n_pos=100,burnin=100,thin = 1,congress = F, hyperparams=list(a = 1, b = 1/10,mu =0, ccc_a = 1, ccc_b=25, kappa_a = 1, omega_sd=0.1, kappa_sd=0.5,
                                                                    i_epi_lower = 0.01, i_epi_upper = 0.08, j_epi_lower = 0.01 ,j_epi_upper = 0.105,
                                                                    i_leap = 10, j_leap = 10,skip = 50, jitter = T, WAIC_group = T),
                                                                    initial_values=NULL,core=10,cluster_seed=8888)
# hyperparams=list(a = 1, b = 1/10, ccc_a = 1, ccc_b=25, kappa_a = 1, omega_sd=0.1, kappa_sd=0.5,
#                  i_epi_lower = 0.005, i_epi_upper = 0.04, j_epi_lower = 0.01 ,j_epi_upper = 0.105,
#                  i_leap = 10, j_leap = 10,skip = 50, jitter = T, WAIC_group = T),
# initial_values=NULL,core=10,cluster_seed=8888)
############



ymat = as.matrix(out[[1]] %>% select(-bioname,-name_district,-icpsr))
pol_info = out[[1]] %>% select(bioname,name_district,icpsr)
dem = grep("\\(D",pol_info$name_district)
gop = grep("\\(R",pol_info$name_district)
ind = grep("\\(I",pol_info$name_district)
rm(out)
pol_info[427,] ## still need to figure out the name district, can manually change it in the end

# init_fun <- function(...) list(mu=4.5, sigma2=0.05)
# 
# library(rstan)
# fit = stan(file = './main_script/Euclidean_1d_stan.stan', data = list(I=nrow(ymat),J = ncol(ymat), y = ymat),seed=8888,control=list(stepsize_jitter = 1,adapt_delta=0.99),iter = 2000,warmup=1000)
#################


congress_member = read_csv('./data/HSall_members.csv',guess_max = 20000)%>% filter(congress == 116)



nomi_nokken = congress_member %>% select(party_code,icpsr,bioname,state_abbrev,district_code,nominate_dim1,nokken_poole_dim1) %>% 
  ## proper formatting of names
  mutate(bioname = str_c(str_sub(bioname,1,1),str_sub(tolower(str_extract(bioname,'[^,]+')),2),str_extract(bioname,', [^,]+'))) %>% 
  ## party code, Democrates: 100, Republican: 200, Independent: 328
  mutate(party = case_when(party_code == 100~'D',party_code == 200~'R',party_code==328~'I',TRUE~'NA')) %>% 
  mutate(name_district  = paste0(bioname," (",party,' ',state_abbrev,"-",district_code,")")) %>% 
  select(icpsr,name_district,bioname,nominate_dim1,nokken_poole_dim1) 





col_stream = rep('red',nrow(pol_info))
col_stream[dem] = 'blue'
col_stream[ind] = 'green'

beta_i_pos = apply(master$beta_i,1,mean)
circ_var = apply(master$beta_i,2,circular::var.circular)
mean(circ_var)

merged_idea_points = pol_info %>% mutate(circular = beta_i_pos) %>% inner_join(nomi_nokken,by=c("name_district"))

x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
plot(round(rank(merged_idea_points$nominate_dim1),digits=0),round(rank(beta_i_pos),digits=0),xlim=c(-5,nrow(pol_info)),ylim=c(-5,nrow(pol_info)),
     col=col_stream, xlab='Euclidean Nominate Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,
     main='',cex.main=2,mgp=c(6,2,0),lwd=2)
abline(1:500,1:500,col='gray',lwd=2)
identify(round(rank(merged_idea_points$nominate_dim1),digits=0),round(rank(beta_i_pos),digits=0),label = pol_info$bioname, cex = 1.5)


x11(height=15,width=15)
par(mar = c(10, 10, 4, 4))
plot(round(rank(merged_idea_points$nokken_poole_dim1),digits=0),round(rank(beta_i_pos),digits=0),xlim=c(-5,nrow(pol_info)),ylim=c(-5,nrow(pol_info)),
     col=col_stream, xlab='Euclidean Nokken Latent Space',ylab='Circular Latent Space',cex.axis=2,cex.lab=2,
     main='',cex.main=2,mgp=c(6,2,0),lwd=2)
abline(1:500,1:500,col='gray',lwd=2)
identify(round(rank(merged_idea_points$nokken_poole_dim1),digits=0),round(rank(beta_i_pos),digits=0),label = pol_info$bioname, cex = 1.5)


