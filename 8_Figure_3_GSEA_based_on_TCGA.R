library(tidyverse)
library(clusterProfiler)
library(enrichplot)
set.seed(1996)

#loading data and custom function for correlation calculation
load("0_data/TCGA_data_extracted.RData")
ensemble_symbol <- readRDS("0_data/list_of_protein_coding_genes_symbol_and_ensembl.RDS")
source('1_functions/column_wise_correlation.R')

# normalization of data, and corelation between genes
expression_tcga <- expression_tcga %>% 
  mutate_at(.vars = 2:ncol(expression_tcga), .funs = ~log2(.+1))

# correlation between COL7A1 and all other available genes must be calculated
#since GSEA will be based on correlation vector
correlation <- column_wise_correlation(set1 = expression_tcga[,-1], 
                                       set2 = "gene", 
                                       gene = "ENSG00000114270.17",
                                       normalization = "none", 
                                       return = "data")

correlation<-correlation %>%  left_join(ensemble_symbol, by = c("set1_genes"="gene_id"))


# performing GSEA based on correlation requires ranked vector
gene_list <- correlation$correlation
names(gene_list) <- correlation$gene_name
gene_list<-sort(gene_list,decreasing = T)


#script used to perform GSEA on few datasets
dbs <- list.files(path = "0_data/MsigDB/", pattern = "[.]gmt$", 
                  recursive = TRUE, full.names = T)

dbs_names <- c("Gene_position","KeggPathways","Transcription_factor",
               "GO_CC","Halmark")

for(d in 1:length(dbs_names)){
  gene_sets <- read.gmt(gmtfile = dbs[d])
  
  
  if(!dir.exists(paths = paste0("CSV/GSEA_MsigDB_TCGA/",dbs_names[d]))){
    dir.create(path = paste0("CSV/GSEA_MsigDB_TCGA/",dbs_names[d]), recursive = T)
  }
  
  if(!dir.exists(paths = paste0("RDS/GSEA_MsigDB_TCGA/",dbs_names[d]))){
    dir.create(path = paste0("RDS/GSEA_MsigDB_TCGA/",dbs_names[d]), recursive = T)
  }
  
  #part needed to transfer gene names

  GSEA_results <- GSEA(geneList = gene_list,
                       TERM2GENE = gene_sets,
                       pvalueCutoff = 1,
                       verbose = T,
                       seed = T,
                       maxGSSize = 2000,
                       nPermSimple = 10000)
  saveRDS(object = GSEA_results, file =  paste0("RDS/GSEA_MsigDB_TCGA/",dbs_names[d],"/", str_sub(str_split_fixed(dbs[d], "/", Inf)[,3], 1, -5), ".RDS")) 
  
  GSEA_results_tibble <- as_tibble(GSEA_results) %>% 
    dplyr::select(Description, ID, NES, p.adjust, everything())
  
  if(dim(GSEA_results_tibble)[1]!=0){
    write_csv(x = GSEA_results_tibble, file = paste0("CSV/GSEA_MsigDB_TCGA/",dbs_names[d],"/", str_sub(str_split_fixed(dbs[d], "/", Inf)[,3], 1, -5), ".csv")) 
  }
}

