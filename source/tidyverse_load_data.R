## function to aggregate legislators' (who switch parties or change district) votes into a single row 
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
get_rollcall_data = function(house,h_s,hn,threshold = 0.4){
  
  ## data path
  PATH = paste0("./data/house/",h_s,hn,"/")
  
  ## loading meta data
  congress_all = read_csv(paste0(PATH,h_s,hn,'_votes.csv'),guess_max = 20000,col_types = cols(prob = col_character()))
  congress_member = read_csv(paste0(PATH,h_s,hn,'_members.csv'),guess_max = 20000)
  congress_bills = read_csv(paste0(PATH,h_s,hn,'_rollcalls.csv'),guess_max = 20000,
                            col_types = cols(session = col_character(),clerk_rollnumber = col_character(),
                                             vote_result = col_character(),vote_question = col_character(),vote_desc = col_character()))
  
  ## create new column name_district by combining name, party affiliation and district
  icpsr_name_district = congress_member %>% select(party_code,icpsr,bioname,state_abbrev,district_code) %>% 
    ## proper formatting of names
    mutate(bioname = str_c(str_sub(bioname,1,1),str_sub(tolower(str_extract(bioname,'[^,]+')),2),str_extract(bioname,', [^,]+'))) %>% 
    ## party code, Democrates: 100, Republican: 200, Independent: 328
    mutate(party = case_when(party_code == 100~'D',party_code == 200~'R',party_code==328~'I',TRUE~'NA')) %>% 
    mutate(name_district  = paste0(bioname," (",party,' ',state_abbrev,"-",district_code,")")) %>% 
    select(icpsr,name_district,bioname) 
  
  ## create data matrix, cast code of 1,2,3 stands for voting "Yea" while 4,5,6 stands for voting "Nah"
  vote_data = congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% 
    mutate(cast_code = case_when(cast_code %in% c(1:3) ~ 1, cast_code %in% c(4:6) ~0,TRUE~NA_real_)) %>% 
    pivot_wider(names_from = rollnumber, values_from = cast_code) %>% select(-congress)
  
  ## joining icpsr with name_district
  full_data = vote_data %>% left_join(icpsr_name_district ,by = "icpsr") %>% 
    relocate(bioname,name_district)
  nrow_full_data = nrow(full_data)
  ## check whether all bills column are valid
  check_bill = full_data %>% select(-bioname,-name_district,-icpsr) %>% colnames() %>% as.numeric() %>% is.na() %>% sum()
  if(check_bill != 0){
    print('Check the column names of the data')
  }
  ##clean memory
  rm(list = c('icpsr_name_district','vote_data','congress_all','congress_member','congress_bills'))
  
  ## find people who switch party or change voting district (have multiple entries/icpsr)
  dup_name_freq = full_data %>% dplyr::count(bioname) %>%  filter(n>1)  %>% select(-n)
  
  ## find legislator who actually change party of voting district within a single congress period
  impossible = full_data %>% dplyr::count(bioname) %>%  filter(n>2) 
  
  if(nrow(impossible)>0){
    print('this legislator is trolling')
    print(impossible)
  }
  
  ## combine these flip floppers' votes into a single row
  if(length(dup_name_freq$bioname)>0){
    full_data_init = full_data %>% filter(!bioname %in% dup_name_freq$bioname)
    for(dup_name in dup_name_freq$bioname){
      full_data_init  = bind_rows(full_data_init ,as_tibble(apply(dup_name_freq %>% inner_join(full_data, by = "bioname") %>% 
                                                               filter(bioname == dup_name) %>% rowwise, 2, aggregate_votes))[1,])
      
    }
    ## check the rows match after combining flip floppers
    check_unique = nrow(full_data_init) == nrow(full_data) - nrow(dup_name_freq)
    if(!check_unique){
      print('check flip floppers')
    }
  }else{
    full_data_init = full_data
    rm(full_data)
  }
  ## number of bills 
  n_bills = full_data_init %>% select(-bioname,-name_district,-icpsr) %>% ncol()
  
  ## check proportion of bills missed in the congress
  missing_prop = apply(full_data_init %>% select(-bioname,-name_district,-icpsr),1,function(x) sum(is.na(x)==T)/n_bills)
  
  ## legislator who missed more than 40% the votes are excluded from the data
  missing_40_percent = full_data_init %>% select(bioname) %>% mutate(missing_prop = missing_prop) %>% filter(missing_prop>0.4)
  
  ## missing more than 40% votes and flip flops!
  troll = missing_40_percent %>% inner_join(dup_name_freq,by='bioname')
  
  ## final data
  full_data_filtered = full_data_init %>% filter(missing_prop <= 0.4)
  rm(full_data_init)
  ## final check preprocessing
  check_row = nrow_full_data == nrow(missing_40_percent) + nrow(full_data_filtered) + nrow(dup_name_freq) - nrow(troll)
  if(!check_row){
    print('Check flip floppers')
  }
  
  return(list(full_data_filtered,dup_name_freq,missing_40_percent))
}

