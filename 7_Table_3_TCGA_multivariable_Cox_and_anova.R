# cox model performed on all expressed "ENSG00000114270.17"s
library(tidyverse)
library(survival)
library("survminer")
library("rms")
load("0_data/TCGA_data_extracted.RData")

# Cox model will avoid using patients for which we have NA for categorical variable
# to make sure that in all cox models we have same patients, we will include only patients
# for which there are no NA values

cox_col7A1 <- expression_tcga %>% 
  select(patient, ENSG00000114270.17) %>% 
  mutate(COL7A1=log2(ENSG00000114270.17+1)) %>% 
  left_join(clinical_tcga) %>% 
  mutate(Stage = factor(Stage, levels = c("Stage I","Stage II", "Stage III","Stage IV"))) %>%
  mutate(Histologic_grade = as.character(Histologic_grade)) %>% 
  mutate(Histologic_grade = ifelse(Histologic_grade %in% c("G1","G2"), "G1/2", Histologic_grade)) %>% 
  mutate(Histologic_grade = factor(Histologic_grade, levels = c("G1/2","G3","G4"))) %>% 
  filter(!(is.na(Histologic_grade) | is.na(Stage)))

#### univariate cox
cox_stage <- coxph(Surv(time = survival_time, event = vital_status)~Stage, data = cox_col7A1)
cox_grade <- coxph(Surv(time = survival_time, event = vital_status)~Histologic_grade ,data = cox_col7A1)
cox_age <- coxph(Surv(time = survival_time, event = vital_status)~Age, data = cox_col7A1)
cox_COL7A1 <- coxph(Surv(time = survival_time, event = vital_status)~COL7A1, data = cox_col7A1)

#summary

summary(cox_stage)
summary(cox_stage)[["sctest"]][["pvalue"]]

summary(cox_grade)
summary(cox_grade)[["sctest"]][["pvalue"]]

summary(cox_age)
summary(cox_age)[["sctest"]][["pvalue"]]

summary(cox_COL7A1)
summary(cox_COL7A1)[["sctest"]][["pvalue"]]


#multivariate cox
cox_clinical <- coxph(Surv(time = survival_time, event = vital_status)~Stage+Histologic_grade+Age, data = cox_col7A1)

cox_multivariable <- coxph(Surv(time = survival_time, event = vital_status)~Stage+Histologic_grade+Age+COL7A1, data = cox_col7A1)
cox_multivariable
summary(cox_multivariable)
anova(cox_multivariable)



  
  