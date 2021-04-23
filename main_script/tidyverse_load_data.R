library(tidyverse)

aggregate_votes = function(x){
  if(sum(is.na(x)) == 2){
    x_mod = NA
  }else{
    x_mod = x[is.na(x) == F]
    if(suppressWarnings(sum(is.na(as.numeric(x_mod)))!=2)){
      x_mod = as.numeric(x_mod)
    }
  }
  return(x_mod)
}

house = T
h_s = ifelse(house==T,'H','S')
hn = 116


get_rollcall_data = function(house,h_s,hn,threshold = 0.4){
  PATH = paste0("./data/house/",h_s,hn,"/")
  
  ## loading meta data
  congress_all = read_csv(paste0(PATH,h_s,hn,'_votes.csv'),guess_max = 20000,col_types = cols(prob = col_character()))
  congress_member = read_csv(paste0(PATH,h_s,hn,'_members.csv'),guess_max = 20000)
  congress_bills = read_csv(paste0(PATH,h_s,hn,'_rollcalls.csv'),guess_max = 20000,
                            col_types = cols(session = col_character(),clerk_rollnumber = col_character(),
                                             vote_result = col_character(),vote_question = col_character(),vote_desc = col_character()))
  
  ## create new column name_district by combining name, party affiliation and district
  ## party code, Democrates: 100, Republican: 200, Independent: 328
  icpsr_name_district = congress_member %>% select(party_code,icpsr,bioname,state_abbrev,district_code) %>% 
    mutate(party = case_when(party_code == 100~'D',party_code == 200~'R',party_code==328~'I',TRUE~'NA')) %>% 
    mutate(name_district  = paste0(bioname," (",party,' ',state_abbrev,"-",district_code,")")) %>% 
    select(icpsr,name_district,bioname) 
  
  ## create data matrix, cast code of 1,2,3 stands for voting Yes while 4,5,6 for voting No 
  vote_data = congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% 
    mutate(cast_code = case_when(cast_code %in% c(1:3) ~ 1, cast_code %in% c(4:6) ~0,TRUE~NA_real_)) %>% 
    pivot_wider(names_from = rollnumber, values_from = cast_code) %>% select(-congress)
  
  
  ## joining icpsr with name_district
  full_data = vote_data %>% left_join(icpsr_name_district ,by = "icpsr") %>% 
    relocate(bioname,name_district)
  
  ##clean memory
  rm(list = c('icpsr_name_district','vote_data','congress_all','congress_member','congress_bills'))
  
  ## find people who switch party or change voting district (have multiple entries/icpsr)
  dup_name_freq = full_data %>% count(bioname) %>%  filter(n>1)  %>% select(-n)
  
  if(length(dup_name_freq$bioname)>0){
    full_data_unique = full_data %>% filter(!bioname %in% dup_name_freq$bioname)
    for(dup_name in dup_name_freq$bioname){
      full_data = bind_rows(full_data_unique,as_tibble(apply(dup_name_freq %>% inner_join(full_data, by = "bioname") %>% 
                                                               filter(bioname == dup_name) %>% rowwise, 2, aggregate_votes))[1,])
      
    }
  }
  return(list(full_data,dup_name_freq))
}

