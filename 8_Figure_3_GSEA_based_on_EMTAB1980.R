library(tidyverse)
library(clusterProfiler)
library(enrichplot)
set.seed(1996)

#loading data and custom function for correlation calculation
load("0_data/EMTAB1980_data_extracted.RData")
source('1_functions/column_wise_correlation.R')

# correlation between COL7A1 and all other available genes must be calculated
#since GSEA will be based on correlation vector
correlation <- column_wise_correlation(set1 = EMTAB1980_expression[,-1], 
                                       set2 = "gene", 
                                       gene = "COL7A1",
                                       normalization = "none", 
                                       return = "data")

# performing GSEA based on correlation requires ranked vector
gene_list <- correlation$correlation
names(gene_list) <- correlation$set1_genes
gene_list<-sort(gene_list,decreasing = T)


#script used to perform GSEA on few datasets
dbs <- list.files(path = "0_data/MsigDB/", pattern = "[.]gmt$", 
                  recursive = TRUE, full.names = T)

dbs_names <- c("Gene_position","KeggPathways","Transcription_factor",
               "GO_CC","Halmark")

for(d in 1:length(dbs_names)){
  gene_sets <- read.gmt(gmtfile = dbs[d])
  
  
  if(!dir.exists(paths = paste0("CSV/GSEA_MsigDB_EMTAB1980/",dbs_names[d]))){
    dir.create(path = paste0("CSV/GSEA_MsigDB_EMTAB1980/",dbs_names[d]), recursive = T)
  }
  
  if(!dir.exists(paths = paste0("RDS/GSEA_MsigDB_EMTAB1980/",dbs_names[d]))){
    dir.create(path = paste0("RDS/GSEA_MsigDB_EMTAB1980/",dbs_names[d]), recursive = T)
  }
  
  #part needed to transfer gene names
  GSEA_results <- GSEA(geneList = gene_list,
                       TERM2GENE = gene_sets,
                       pvalueCutoff = 1,
                       verbose = T,
                       seed = T,
                       maxGSSize = 2000,
                       nPermSimple = 10000)
  
  saveRDS(object = GSEA_results, file =  paste0("RDS/GSEA_MsigDB_EMTAB1980/",dbs_names[d],"/", str_sub(str_split_fixed(dbs[d], "/", Inf)[,3], 1, -5), ".RDS")) 
  
  GSEA_results_tibble <- as_tibble(GSEA_results) %>% 
    dplyr::select(Description, ID, NES, p.adjust, everything())
  
  if(dim(GSEA_results_tibble)[1]!=0){
    write_csv(x = GSEA_results_tibble, file = paste0("CSV/GSEA_MsigDB_EMTAB1980/",dbs_names[d],"/", str_sub(str_split_fixed(dbs[d], "/", Inf)[,3], 1, -5), ".csv")) 
  }
}



