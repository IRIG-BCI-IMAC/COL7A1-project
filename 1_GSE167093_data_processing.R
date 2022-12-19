#### wozniak data extraction and normalization 
library(tidyverse)
library(readxl)
 # GSE167093 dataset

GSE167093_expression <- read_excel("0_data/GSE167093/GSE167093_Tumor_only_norm_matrix.xlsx")
GSE167093_clinical <- read_delim(file='0_data/GSE167093/wozniak reworked.tsv', comment = "#")
GSE167093_array_explain <- read_delim(file='0_data/GSE167093/GPL10558-50081.txt', comment = "#")

#permuting the data to get the needed format and average values 
colnames(GSE167093_expression)[1] <- "ID"

GSE167093_array_explain<- GSE167093_array_explain %>% 
  select(ID, Symbol) %>% 
  filter(!is.na(Symbol)) %>% 
  filter(Symbol != "permuted_negative")


GSE167093_expression <- GSE167093_expression %>% 
  left_join(GSE167093_array_explain) %>% 
  select(-ID) %>% 
  select(Symbol, everything()) %>% 
  group_by(Symbol) %>% 
  summarise_all(.funs = mean) %>% 
  ungroup()%>% 
  pivot_longer(cols = 2:605, values_to= "expression", names_to="patient") %>% 
  pivot_wider(names_from = Symbol, values_from=expression)

############# part where I load TCGA data and make sure same genes are involved

save(list =c("GSE167093_expression","GSE167093_clinical"), file = "0_data/GSE167093_data_extracted.RData" )
