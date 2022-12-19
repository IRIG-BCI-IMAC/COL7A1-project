# cox model performed on all expressed "COL7A1"s
library(tidyverse)
library(survival)
library(survminer)

load("0_data/TCGA_data_extracted.RData")
source(file = "1_functions/kaplan_meier_curve.R")

if(!dir.exists(paths = "Figure_2/")){
  dir.create(path = "Figure_2/", recursive = T)
}


expression_tcga <-expression_tcga %>% 
  select(patient, ENSG00000114270.17 ) %>% 
  mutate(COL7A1=log2(ENSG00000114270.17+1))

# Using cutoff fucntion to determine best spot to cut dataset in 2 groups
cutoff_df <- tibble(time=clinical_tcga$survival_time,
                    event=clinical_tcga$vital_status,
                    COL7A1=expression_tcga$COL7A1)
cutpoint <- surv_cutpoint(cutoff_df, variables = "COL7A1")
plot(cutpoint, "COL7A1", palette = "npg")

## ## Using function to remake Kaplan-Meier curves for TCGA
## stage I
stage_I <- clinical_tcga %>% 
  filter(Stage=="Stage I")
stage_I<-expression_tcga %>% 
  select(patient, "COL7A1") %>% 
  semi_join(stage_I)


KM_TCGA_COL7A1_stage_1 <- KM_curve(mRNA_counts = stage_I, clinical_data = clinical_tcga, gene = "COL7A1",
                                   grouping = "set range", censoring = T, use_outliers = T,
                                   group_names=c("Low","High"),
                                   ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,5), 
                                   histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))

svg("Figure_2/Stage_I.svg")
KM_TCGA_COL7A1_stage_1[[2]]
dev.off()
## stage II

stage_II <- clinical_tcga %>% 
  filter(Stage=="Stage II")
stage_II<-expression_tcga %>% 
  select(patient, "COL7A1") %>% 
  semi_join(stage_II)

KM_TCGA_COL7A1_stage_2 <- KM_curve(mRNA_counts = stage_II, clinical_data = clinical_tcga, gene = "COL7A1",
                                   grouping = "set range", censoring = T, use_outliers = T,
                                   group_names=c("Low","High"),
                                   ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,5), 
                                   histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))

svg("Figure_2/Stage_II.svg")
KM_TCGA_COL7A1_stage_2[[2]]
dev.off()
## stage III

stage_III <- clinical_tcga %>% 
  filter(Stage=="Stage III")
stage_III<-expression_tcga %>% 
  select(patient, "COL7A1") %>% 
  semi_join(stage_III)

KM_TCGA_COL7A1_stage_3 <- KM_curve(mRNA_counts = stage_III, clinical_data = clinical_tcga, gene = "COL7A1",
                                   grouping = "set range", censoring = T, use_outliers = T,
                                   group_names=c("Low","High"),
                                   ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,5), 
                                   histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))


svg("Figure_2/Stage_III.svg")
KM_TCGA_COL7A1_stage_3[[2]]
dev.off()
## stage IV

stage_IV <- clinical_tcga %>% 
  filter(Stage=="Stage IV")
stage_IV<-expression_tcga %>% 
  select(patient, "COL7A1") %>% 
  semi_join(stage_IV)

KM_TCGA_COL7A1_stage_4 <- KM_curve(mRNA_counts = stage_IV, clinical_data = clinical_tcga, gene = "COL7A1",
                                   grouping = "set range", censoring = T, use_outliers = T,
                                   group_names=c("Low","High"),
                                   ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,5), 
                                   histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))

svg("Figure_2/Stage_IV.svg")
KM_TCGA_COL7A1_stage_4[[2]]
dev.off()
# histogram
hist_data<-clinical_tcga %>% filter(Stage!="")
hist_data <- expression_tcga %>% 
  semi_join(hist_data) 

histogram <- KM_curve(mRNA_counts = hist_data, clinical_data = clinical_tcga, gene = "COL7A1",
                      grouping = "set range", censoring = T, use_outliers = T,
                      group_names=c("Low","High"), histogram_bins=20,
                      ranges=cutpoint$cutpoint$cutpoint,facet_by = "Stage",
                      histogram_limits_x = c(0,7.5), 
                      histogram_limits_y = c(0,1.7),  survival_limit_x=c(0,10))


svg("Figure_2/Stage_split_histogram.pdf", width = 28)
histogram[[1]]
dev.off()


