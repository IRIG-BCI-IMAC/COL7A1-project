# creating the KM curve and general COL7A1 distribution
library(tidyverse)
library(survival)

load("0_data/EMTAB1980_data_extracted.RData")
source(file = "1_functions/kaplan_meier_curve.R")


if(!dir.exists(paths = "Figure_1/")){
  dir.create(path = "Figure_1/", recursive = T)
}



cutoff_df <- tibble(time=EMTAB1980_clinical$survival_time,
                    event=EMTAB1980_clinical$vital_status,
                    COL7A1=EMTAB1980_expression$COL7A1)
cutpoint <- surv_cutpoint(cutoff_df, variables = "COL7A1",minprop = 0.1)

plot(cutpoint, "COL7A1", palette = "npg")

svg(file = "Figure_1/EMTAB1980_KM_2_levels.svg")

KM_curve(mRNA_counts = EMTAB1980_expression,
         clinical_data = EMTAB1980_clinical,
         grouping="set range",ranges = cutpoint$cutpoint$cutpoint,
         gene = "COL7A1",histogram_bins = 10,
         use_outliers = F,
         group_names =factor(c("1:Low","2:High")))[[2]]
dev.off()



