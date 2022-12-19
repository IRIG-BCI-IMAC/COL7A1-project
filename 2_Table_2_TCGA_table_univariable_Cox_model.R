# cox model performed on all expressed genes
library(tidyverse)
library(survival)


# Data loading ------------------------------------------------------------
load("0_data/TCGA_data_extracted.RData")

#in case I need symbols for final table
ensembl_to_symbol <- read_tsv(file = "GDCdata/TCGA-KIRC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/000d81dd-9ba4-4852-9090-2bf22f6483f0/09755ce8-ed89-411a-a42f-b3edc4e41eeb.rna_seq.augmented_star_gene_counts.tsv", comment = "N_", skip = 1)[,c(1,2,3)]


survival_data <- as.matrix(clinical_tcga[,c(2,7)], )
rownames(survival_data) <- clinical_tcga$patient


# Data transformation -----------------------------------------------------

gene_data <- expression_tcga %>% 
  mutate_at(.vars = 2:ncol(expression_tcga), .funs = ~log2(.+1))
barcodes <- gene_data$patient

gene_data<- as.matrix(gene_data[,-1])
rownames(gene_data) <- barcodes
rm(expression_tcga)


# Cox model  ------------------------------------------------------------

#Runing the cox model on all genes found in dataset, as well as extracting 
#the data about individual genes

table_cox=NULL

surv_data <- Surv(survival_data[,2], survival_data[,1])

for( i in 2:ncol(gene_data)){
  #### Cox model and basic statistics ####
  gene_vect <- gene_data[,i]
  
  median <- median(gene_vect)
  MAD<- mad(gene_vect)
  quan<- quantile(x = gene_vect, probs = c(0.2,0.8))
  Q20<-unname(quan[1])
  Q80 <- unname(quan[2])
  R_80_20=Q80-Q20
  cox_data<- summary(coxph(formula =surv_data~gene_vect))
  cox_gene=tibble(ENSEMBL_id=colnames(gene_data)[i],
                  p.value=cox_data[["sctest"]][["pvalue"]],
                  LogTest=cox_data[["sctest"]][["test"]],
                  beta=cox_data$coefficients[1],
                  HR=cox_data$coefficients[2],
                  median=median,
                  MAD=MAD,
                  Q20=Q20, 
                  Q80=Q80, 
                  R_80_20=R_80_20)
  table_cox <- rbind(table_cox,cox_gene)
}

# to filter out genes that are almost non-expressed

table_cox<-table_cox %>% 
  left_join(ensembl_to_symbol, by=c("ENSEMBL_id" = "gene_id")) %>% 
  filter(!median==0 & !R_80_20==0) %>% 
  #select(ENSEMBL_id, gene_name, everything()) %>% 
  arrange(p.value) %>% 
  mutate(rank=row_number())

if(!dir.exists(paths = "CSV/1_Cox_model_table/")){
  dir.create(path = "CSV/1_Cox_model_table/", recursive = T)
}

write_csv(x = table_cox, file = "CSV/1_Cox_model_table/model_on_all_genes.csv")




