# creating the KM curve and general COL7A1 distribution
library(tidyverse)
library(survival)

load("0_data/TCGA_data_extracted.RData")
source(file = "1_functions/kaplan_meier_curve.R")


if(!dir.exists(paths = "Figure_1/")){
  dir.create(path = "Figure_1/", recursive = T)
}

KM_data_expression <- expression_tcga %>% 
  select(patient, ENSG00000114270.17) %>% 
  mutate(COL7A1 = log2(ENSG00000114270.17+1))

svg(file = "Figure_1/TCGA_KM_4_levels.svg")
KM_curve(mRNA_counts = KM_data_expression,
         clinical_data = clinical_tcga,
         group_names = factor(c("1","2","3","4")),
         gene = "COL7A1",
         use_outliers = F)[[2]]


dev.off()

svg(file = "Figure_1/TCGA_KM_4_levels_hist.svg")
KM_curve(mRNA_counts = KM_data_expression,
         clinical_data = clinical_tcga,
         group_names = factor(c("1","2","3","4")),
         gene = "COL7A1",histogram_bins = 2,
         use_outliers = F)[[1]]


dev.off()


cutoff_df <- tibble(time=clinical_tcga$survival_time,
                    event=clinical_tcga$vital_status,
                    COL7A1=KM_data_expression$COL7A1)
cutpoint <- surv_cutpoint(cutoff_df, variables = "COL7A1")

plot(cutpoint, "COL7A1", palette = "npg")

svg(file = "Figure_1/TCGA_KM_2_levels.svg")
KM_curve(mRNA_counts = KM_data_expression,
         clinical_data = clinical_tcga,
         grouping="set range",ranges = cutpoint$cutpoint$cutpoint,
         gene = "COL7A1",
         use_outliers = F, histogram_bins = 4,
         group_names =factor(c("1:Low","2:High")))[2]
dev.off()


