#### EMTAB1980 data extraction and normalization 
library(tidyverse)

EMTAB1980_expression <- read_delim(file='0_data/E_MTAB_1980/ccRCC_exp_log_quantile_normalized.txt')
EMTAB1980_clinical <- read_delim(file='0_data/E_MTAB_1980/ccRCC_clinical.csv')
EMTAB1980_array_explain <- read_delim(file='0_data/E_MTAB_1980/A-MEXP-2183_comments.txt')


EMTAB1980_array_explain<- EMTAB1980_array_explain %>% 
  distinct(`Comment[GeneName]`,`Comment[SystematicName]`) %>% 
  rename('Symbol'="Comment[GeneName]",'SystematicName'="Comment[SystematicName]")


EMTAB1980_expression <- EMTAB1980_expression %>% 
  select(-1,-2) %>%
  left_join(EMTAB1980_array_explain) %>% 
  select(-SystematicName) %>% 
  select(Symbol, everything()) %>% 
  group_by(Symbol) %>% 
  summarise_all(.funs = mean) %>% 
  ungroup()%>% 
  pivot_longer(cols = 2:102, values_to= "expression", names_to="patient") %>% 
  pivot_wider(names_from = Symbol, values_from=expression)

# write_csv(x = sat_average, file = "csv_output/EMTAB1980_expression.csv")
#### EMTAB1980 clinical data 


EMTAB1980_clinical <- EMTAB1980_clinical %>% 
  select(`sample ID`,Sex, Age,`Stage at diagnosis`,`Fuhrman grade`,outcome,`observation period (month)`) %>% 
  rename("stage"="Stage at diagnosis") %>% 
  mutate(t_grade=str_sub(stage, start=2, end=3)) %>% 
  mutate(n_grade=str_sub(stage, start=-4, end=-3)) %>% 
  mutate(m_grade=str_sub(stage, start=-2, end=-1)) %>% 
  mutate(Stage= case_when(t_grade=="T1" & n_grade=="N0" & m_grade=="M0" ~ "stage I",
                          t_grade== "T2" & n_grade=="N0" & m_grade=="M0" ~ "stage II",
                          t_grade== "T3" & n_grade=="N0" & m_grade=="M0" ~ "stage II",
                          n_grade=="N1" & m_grade=="M0" ~ "stage III",
                          n_grade=="N2" & m_grade=="M0" ~ "stage III",
                          t_grade== "T4"  & m_grade=="M0" ~ "stage III",
                          t_grade== "T1" & n_grade=="N0" & m_grade=="M0" ~ "stage I",
                          m_grade=="M1" ~ "stage IV")) %>% 
  rename("patient"="sample ID") %>% 
  rename("Histologic_grade"="Fuhrman grade") %>% 
  rename("survival_time"="observation period (month)") %>% 
  rename("vital_status"="outcome") %>% 
  mutate(survival_time = survival_time/12) %>% 
  mutate(vital_status = case_when(vital_status=="alive" ~ 0,
                                  vital_status=="dead" ~ 1)) %>% 
  mutate(patient= str_replace(patient, pattern = '\\.', replacement = '-')) %>% 
  filter(patient %in% EMTAB1980_expression$patient)


# saving the data

save(list =c("EMTAB1980_expression","EMTAB1980_clinical"), file = "0_data/EMTAB1980_data_extracted.RData" )