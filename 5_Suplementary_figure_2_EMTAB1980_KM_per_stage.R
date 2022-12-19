# cox model performed on all expressed "ENSG00000114270.17"s
library(tidyverse)
library(survival)
library(survminer)

load("0_data/EMTAB1980_data_extracted.RData")
source(file = "1_functions/kaplan_meier_curve.R")

if(!dir.exists(paths = "Suplementary_figure_2/")){
  dir.create(path = "Suplementary_figure_2/", recursive = T)
}

# Using cutoff fucntion to determine best spot to cut dataset in 2 groups
cutoff_df <- tibble(time=EMTAB1980_clinical$survival_time,
                    event=EMTAB1980_clinical$vital_status,
                    COL7A1=EMTAB1980_expression$COL7A1)
cutpoint <- surv_cutpoint(cutoff_df, variables = "COL7A1")
plot(cutpoint, "COL7A1", palette = "npg")

## ## Using function to remake Kaplan-Meier curves for TCGA
## stage I
stage_I_II <- EMTAB1980_clinical %>% 
  filter(Stage=="stage I" |Stage=="stage II")
stage_I_II<-EMTAB1980_expression %>% 
  select(patient, "COL7A1") %>% 
  filter(patient %in% stage_I_II$patient)


KM_TCGA_COL7A1_stage_1_2 <- KM_curve(mRNA_counts = stage_I_II, 
                                     clinical_data = EMTAB1980_clinical, 
                                     gene = "COL7A1",
                                     grouping = "set range", censoring = T, 
                                     group_names=c("Low","High"),
                                     ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,10), 
                                     histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))
## stage II

stage_III_IV <- EMTAB1980_clinical %>% 
  filter(Stage=="stage III" |Stage=="stage IV")
stage_III_IV<-EMTAB1980_expression %>% 
  select(patient, "COL7A1") %>% 
  filter(patient %in% stage_III_IV$patient)

KM_TCGA_COL7A1_stage_3_4 <- KM_curve(mRNA_counts = stage_III_IV, clinical_data = EMTAB1980_clinical, gene = "COL7A1",
                                     grouping = "set range", censoring = T, use_outliers = T,
                                     group_names=c("Low","High"),
                                     ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,10), 
                                     histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))

# histogram
hist_data<-EMTAB1980_clinical %>% 
  filter(Stage!="") %>% 
  mutate(Stage=ifelse(Stage %in% c("stage I","stage II"),"Stage I/II","Stage III/IV"))
EMTAB1980_expression <- EMTAB1980_expression %>% 
  filter(patient %in% hist_data$patient) 

histogram <- KM_curve(mRNA_counts = EMTAB1980_expression, clinical_data = hist_data,
                      gene = "COL7A1",
                      grouping = "set range", censoring = T, use_outliers = T,
                      group_names=c("Low","High"), histogram_bins=5,
                      ranges=cutpoint$cutpoint$cutpoint,facet_by = "Stage",
                      histogram_limits_x = c(3,10), 
                      histogram_limits_y = c(0,1),  survival_limit_x=c(0,10))


svg("Suplementary_figure_2/Stage_split_1+2_EMTAB1980.svg")
KM_TCGA_COL7A1_stage_1_2[[2]]+ggtitle("Stage I/II")
dev.off()

svg("Suplementary_figure_2/Stage_split_3+4_EMTAB1980.svg")
KM_TCGA_COL7A1_stage_3_4[[2]]+ggtitle("Stage III/IV")
dev.off()

svg("Suplementary_figure_2/Stage_histogram_EMTAB1980.svg", width = 15)
histogram[[1]]
dev.off()


