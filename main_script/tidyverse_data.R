library(tidyverse)

congress_all = read_csv('./data/Hall_votes.csv',guess_max = 20000)
congress_member = read_csv('./data/HSall_members.csv',guess_max = 20000)
congress_bills = read_csv('./data/Hall_rollcalls.csv',guess_max = 20000,col_types = cols(session = col_character(),
 clerk_rollnumber = col_character(),vote_result = col_character(),vote_question = col_character(),vote_desc = col_character()))

H116_bills = read_csv('./data/H116_rollcalls.csv',guess_max =20000)

congress_all %>% glimpse()
congress_member %>% glimpse()
H116_bills %>% glimpse()

H116_bills$bill_number
###############################
congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress %in% c(100:116)) %>%
  count(cast_code)
## only three parties in 100 to 116
congress_member %>% filter(congress %in% c(100:116)) %>%
  count(party_code)

congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress == 116) %>%
  count(cast_code)

h116_pol_info = congress_member %>% filter(congress == 116) %>% select(party_code,icpsr,bioname,state_abbrev,district_code) %>% 
  mutate(party = case_when(party_code == 100~'D',party_code == 200~'R',party_code==328~'I')) %>% 
  mutate(name_district  = paste0(bioname," (",party,' ',state_abbrev,"-",district_code,")")) %>% 
  select(icpsr,name_district,bioname) 

###############################
h116 = congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress == 116) %>% 
  mutate(cast_code = case_when(cast_code %in% c(1:3) ~ 1, cast_code %in% c(4:6) ~0,TRUE~NA_real_)) %>% 
  # inner_join(H116_bills %>% select(rollnumber,))
  pivot_wider(names_from = rollnumber, values_from = cast_code) %>% select(-congress)

apply(h116 %>% select(-icpsr),2,function(x) sum(x,na.rm=T)) == H116_bills$yea_count

apply(h116 %>% select(-icpsr),2,function(x) sum(x,na.rm=T)) == H116_bills$yea_count


## check yea and nah count 
temp = h116 %>% select(-icpsr) %>% summarise(across(everything(),function(x) sum(x,na.rm = T)))
near(as.numeric(unlist(temp)),H116_bills$yea_count)
temp2= h116 %>% select(-icpsr) %>% summarise(across(everything(),function(x) length(which(x==0))))
near(as.numeric(unlist(temp2)),H116_bills$nay_count)


congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress == 116) %>% 
  summarise(count = n_distinct(rollnumber))

h116_all = h116 %>% left_join(h116_pol_info,by = "icpsr")



temp = h116_all %>% select(-bioname,-name_district,-icpsr)
temp2 = apply(temp ,1,function(x) length(which(is.na(x)==T))) / ncol(temp )
temp3 = which(temp2>0.4)

h116_all %>% group_by(bioname) %>%  mutate(n = n()) %>% filter(n>1) %>% 
  select(bioname,name_district)

h116_all %>% count(bioname) %>% filter(n>1)
h116_all[temp3,] %>% select(bioname)

house_list = c(100:116)
dup_list = yes_no_list = vector('list',17)

for(i in 1:17){
  cat("\rProgress: ",i,"/",17)
  house_number = house_list[i]
  
  c1 = congress_member %>% filter(congress == house_number) %>% select(party_code,icpsr,bioname,state_abbrev,district_code) %>% 
    mutate(party = case_when(party_code == 100~'D',party_code == 200~'R',party_code==328~'I')) %>% 
    mutate(name_district  = paste0(bioname," (",party,' ',state_abbrev,"-",district_code,")")) %>% 
    select(icpsr,name_district,bioname) 
  
  c2= congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress == house_number) %>% 
    mutate(cast_code = case_when(cast_code %in% c(1:3) ~ 1, cast_code %in% c(4:6) ~0,TRUE~NA_real_)) %>% 
    pivot_wider(names_from = rollnumber, values_from = cast_code) %>% select(-congress)
  
  c_bills = congress_bills %>% filter(congress == house_number)
    
  temp = c2 %>%select(-icpsr) %>%  summarise(across(everything(),function(x) sum(x,na.rm = T)))
  t1 = sum(near(as.numeric(unlist(temp)),c_bills$yea_count)) == length(temp)
  
  ###H115 anomaly for yes_no_list, after checking, it's trump's votes! 99912 is trump
  # as_tibble(as.numeric(unlist(temp)) - c_bills$yea_count) %>% mutate(rollnumber=c_bills$rollnumber) %>%
  #   filter(value!=0) %>% print(n=Inf) %>% inner_join(  congress_all %>% select(icpsr,cast_code,congress,rollnumber)
  #                                                      %>% filter(congress == house_number),by=c("rollnumber") ) %>%
  #   select(-icpsr,-congress) %>% group_by(rollnumber) %>% count(cast_code) %>%
  #   pivot_wider(names_from = cast_code,values_from = n) %>% print(n = Inf) %>% inner_join(c_bills %>%
  #                                     select(rollnumber,yea_count),by=c("rollnumber")) %>% print(n=Inf) %>%
  # inner_join(c2 %>% filter(icpsr==99912) %>% pivot_longer(cols=-icpsr) %>% rename(rollnumber = name) %>%
  #              mutate(rollnumber = as.numeric(rollnumber) )%>% select(-icpsr),
  #  by=c("rollnumber")) %>% mutate(corrected_yea = yea_count + value) %>% 
  #   mutate(check = `1`==corrected_yea)%>% print(n = Inf)
  #   
    ###
  temp2 = c2 %>%select(-icpsr) %>%  summarise(across(everything(),function(x) length(which(x==0))))
  t2 = sum(near(as.numeric(unlist(temp2)),c_bills$nay_count)) == length(temp)
  yes_no_list[[i]] = c(t1,t2)
  
  c3 = c2 %>% left_join(c1,by = "icpsr")
  dup_list[[i]] = c3 %>% group_by(bioname) %>%  mutate(n = n()) %>% filter(n>1) %>% 
    select(bioname,name_district) %>% arrange(desc(name_district))
  ### work here
  c3 %>% group_by(bioname) %>%  mutate(n = n()) %>% filter(n>1) %>% select(bioname)

  ###########
  c3 %>% group_by(bioname) %>%  mutate(n = n()) %>% filter(n>1) %>% 
    select(bioname,name_district) %>% arrange(desc(name_district))
  ###########
  dup_list[[i]]$name_district = gsub("\\)","",gsub(".*\\(","",dup_list[[i]]$name_district))
  
  rm(c1)
  rm(c2)
  rm(c3)
  rm(c_bills)
  rm(t1)
  rm(t2)
  rm(temp2)
  rm(temp)
}
names(dup_list) = names(yes_no_list) = paste0('H',100:116)


dup_list_pivot = lapply(dup_list,function(x) x %>% group_by(bioname) %>% mutate(row = row_number()) %>% 
  pivot_wider(names_from = bioname,values_from=name_district) %>% select(-row))

yes_no_list




library(wnominate)


vote2 = readKH2(file='F:\\Download2\\h116_votes.ord')
 
vote<-as.matrix(vote2$votes)
####Amash switches party at 430, we combine the record for Amash together####
#   vote[202,430:ncol(vote)] = vote[203,430:ncol(vote)]
#   vote = vote[-203,]
# }else{
#   vote<-as.matrix(vote2$votes)
# }
pol<-rownames(vote)
####0,7,8,9 represents NA###
ind3<-which(vote==0 | vote==7 | vote==8| vote==9)
vote[ind3]<-NA
####1,2,3 represents Yea####
ind1<-which(vote==1 | vote==2 | vote==3)
vote[ind1]<-1
####4,5,6 represents Nah####
ind2<-which(vote==4 | vote==5| vote==6)
vote[ind2]<-0

###tidyverse gives same result
sum(vote,na.rm = T) == h116 %>% select(-icpsr) %>% sum(na.rm=T)



test_df = matrix(rnorm(400*1500),400,1500)

test_tibble = as_tibble(test_df)

x = rnorm(1500)




