#NEW dataset

library(tidyverse)
library(readxl)
library(survminer)
library(survival)
source(file = "1_functions/kaplan_meier_curve.R")


# I will extract OS and event data from fist sheed of suplementary table
#from expression sheet I will extract expression data for cohort CM-025 that was
#treated with NIVOLUMAB


clinical_braun <- read_excel("0_data/Braun_data/41591_2020_839_MOESM2_ESM.xlsx", 
                       sheet = "S1_Clinical_and_Immune_Data", 
                       skip = 1)
expression_braun <- read_excel("0_data/Braun_data/41591_2020_839_MOESM2_ESM.xlsx", 
                         sheet = "S4A_RNA_Expression", skip = 1)


clinical_braun <- clinical_braun %>% 
  filter(RNA_ID != "NA") %>% 
  filter(Tumor_Sample_Primary_or_Metastasis == "PRIMARY") %>% 
  arrange(RNA_ID)%>% 
  rename("survival_time"="OS") %>% 
  mutate(survival_time=survival_time/12) %>% 
  rename("vital_status"="OS_CNSR") %>% 
  mutate(patient=RNA_ID)



ensemble_symbol <- readRDS("0_data/list_of_protein_coding_genes_symbol_and_ensembl.RDS")

expression_braun <- expression_braun %>% 
  filter(gene_name %in% ensemble_symbol$gene_name) %>% 
  select(gene_name, any_of(clinical_braun$RNA_ID)) %>% 
  mutate_at(2:226, .funs = ~.-21) %>% 
  pivot_longer(cols = 2:226, names_to = "patient", values_to = "expression") %>% 
  pivot_wider(names_from = gene_name, values_from = expression) %>% 
  arrange(patient)

save(list = c("clinical_braun", "expression_braun"), file = "0_data/Braun_data_extracted.RData")
