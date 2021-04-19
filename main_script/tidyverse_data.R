library(tidyverse)

congress_all = read_csv('./data/Hall_votes.csv',guess_max = 20000)
congress_member = read_csv('./data/HSall_members.csv',guess_max = 20000)


congress_all %>% glimpse()
congress_member %>% glimpse()


h116 = congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress == 116) %>% 
  pivot_wider(names_from = rollnumber, values_from = cast_code) %>% select(-congress)

congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress == 116) %>%
 count(cast_code)

congress_all %>% select(icpsr,cast_code,congress,rollnumber) %>% filter(congress %in% c(100:116)) %>%
  count(cast_code)
