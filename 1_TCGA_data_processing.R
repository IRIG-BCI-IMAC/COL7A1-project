# Loading packages --------------------------------------------------------

library(TCGAbiolinks)
library(tidyverse)

# Preparing a query for ccRCC data: clinical and RNAseq -------------------
# needed to load and download data. run every time
query <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts")

query_clinical <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML")


# using TCGAbiolinks to download unfiltered data --------------------------
# Data should be downloaded only once. I will leave functions as comments
#GDCdownload(query_clinical)
#GDCdownload(query = query)


# loading clinical data ---------------------------------------------------
clinical <- GDCquery_clinic("TCGA-KIRC","clinical")
clinical_additional <- GDCprepare_clinic(query_clinical,clinical.info = "patient") %>% 
  select(bcr_patient_barcode,neoplasm_histologic_grade)
clinical<- left_join(clinical, clinical_additional)

# loading of RNAseq data --------------------------------------------------
files_data_primary_tumor <- getResults(query) %>%  #file_id and barcodes needed
  filter(sample_type=="Primary Tumor")

#getting the path to file and reading first file to obtaind 
#ensemble ID and gene_name and gene type of genes/rows. This will serve as scafold for automated
# data loading. in For loop I will just add patient by patient expression data.

file_read<-list.files(path = paste0("GDCdata/TCGA-KIRC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/",files_data_primary_tumor[1,1]), full.names = T)
data <- read_tsv(file = file_read, comment = "N_", skip = 1)

primary_tumor_tpm <- read_tsv(file = file_read, comment = "N_", skip = 1)[,c(1,2,3)]
rm(data)
# reading all primary tumor tissue samples, extracting the TPM counts,
#and using the barcode as colname for each patient

for(i in 1:dim(files_data_primary_tumor)[1]){
  file_read<-list.files(path = paste0("GDCdata/TCGA-KIRC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/",files_data_primary_tumor[i,1]), full.names = T)
  data <- read_tsv(file = file_read, comment = "N_", skip = 1)[,c(1,2,7)]
  primary_tumor_tpm <- left_join(primary_tumor_tpm,data, by=c("gene_id", "gene_name"))
  colnames(primary_tumor_tpm)[i+3] <- files_data_primary_tumor[i,3]
}


#save(list = c("primary_tumor_tpm","clinical"), file = "KIRC_tpm_and_clinical.RDS")


# Filtering patient based on Firebrowse, gene type and expression---------------------------

# I have noticed that number of samples is higher than in FireBrowse data used before
# few samples are for the same patient, and are shown as outliers on PCA.
# I will use colnames of FireBrowse KIRC data to select same patients/samples.
# I will filter out nonexpressed genes (TPM < 10 in all samples) and non-protein-coding genes

samples_in_firebrowse <- readRDS("0_data/patients_in_firebrowse.RDS")
primary_tumor_tpm <-primary_tumor_tpm[,c(T,T,T,colnames(primary_tumor_tpm)[c(-1,-2,-3)] %in% samples_in_firebrowse)]

# summ of counts per gene 
expression_selector <- apply(primary_tumor_tpm[,-(1:3)], 1, sum) > 10
gene_type_selector <- primary_tumor_tpm$gene_type=="protein_coding"
primary_tumor_tpm <- primary_tumor_tpm[expression_selector&gene_type_selector,]

primary_tumor_tpm_protein <- primary_tumor_tpm[,-(2:3)] %>% 
  pivot_longer(cols = 2:(ncol(primary_tumor_tpm)-2), values_to= "expression", names_to="patient") %>% 
  pivot_wider(names_from = "gene_id", values_from='expression')

# clearing up barcode so it would relate to patients from clinical data 

primary_tumor_tpm_protein$patient <- str_sub(primary_tumor_tpm_protein$patient, 1, 12)

# clinical data filtering -------------------------------------------------

clinical <- clinical[clinical$bcr_patient_barcode %in% primary_tumor_tpm_protein$patient,] %>% 
  select(bcr_patient_barcode, vital_status,
         days_to_death,
         days_to_last_follow_up, 
         neoplasm_histologic_grade,
         ajcc_pathologic_stage,
         age_at_index,
         gender) %>% 
  mutate(survival_time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death)) %>% 
  mutate(vital_status = ifelse(vital_status=="Dead", 1,0)) %>% 
  select(-days_to_last_follow_up, -days_to_death) %>% 
  mutate(survival_time = survival_time/365.25) %>% 
  rename(c("patient"="bcr_patient_barcode",
           "Histologic_grade"="neoplasm_histologic_grade",
           "Stage"="ajcc_pathologic_stage",
           "Age"="age_at_index",
           "Sex"="gender"))

# Sorting and saving data -------------------------------------------------
expression_tcga <- primary_tumor_tpm_protein %>% arrange(patient)

clinical_tcga <- clinical %>% arrange(patient)

save(list = c("expression_tcga","clinical_tcga"), file = "0_data/TCGA_data_extracted.RData")
