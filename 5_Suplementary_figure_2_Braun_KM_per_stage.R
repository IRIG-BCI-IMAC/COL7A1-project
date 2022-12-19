# creating the KM curve and general COL7A1 distribution
library(tidyverse)
library(survival)

load("0_data/Braun_data_extracted.RData")
source(file = "1_functions/kaplan_meier_curve.R")


if(!dir.exists(paths = "Suplementary_figure_2/")){
  dir.create(path = "Suplementary_figure_2/", recursive = T)
}



cutoff_df <- tibble(time=clinical_braun$survival_time,
                    event=clinical_braun$vital_status,
                    COL7A1=expression_braun$COL7A1)
cutpoint <- surv_cutpoint(cutoff_df, variables = "COL7A1")
cutoff_df <- surv_categorize(cutpoint, variables = "COL7A1")
plot(cutpoint, "COL7A1", palette = "npg")

braun_km <- KM_curve(mRNA_counts = expression_braun,
                     clinical_data = clinical_braun,
                     grouping ="set range",
                     ranges = cutpoint$cutpoint$cutpoint,
                     gene = "COL7A1",
                     use_outliers = F,
                     histogram_bins = 5,
                     group_names = factor(c("1:Low","2:High")),
                     histogram_limits_x = c(5,17), 
                     histogram_limits_y = c(0,1), 
                     survival_limit_x=c(0,10))

svg(file = "Suplementary_figure_2/Stage4_KM_curve_braun.svg")
braun_km[[2]]
dev.off()

svg(file = "Suplementary_figure_2/Stage4_KM_curve_braun_hist.svg")
braun_km[[1]]
dev.off()

# Division of expression according to treatment


braun_km_niv <- KM_curve(mRNA_counts = expression_braun,
                         clinical_data = subset(clinical_braun, Arm=="NIVOLUMAB"),
                         grouping ="set range",
                         ranges = cutpoint$cutpoint$cutpoint,
                         gene = "COL7A1",
                         use_outliers = F,
                         histogram_bins = 5,
                         group_names = factor(c("1:Low","2:High")),
                         histogram_limits_x = c(5,17), 
                         histogram_limits_y = c(0,1), 
                         survival_limit_x=c(0,10))

braun_km_evr<- KM_curve(mRNA_counts = expression_braun,
                         clinical_data = subset(clinical_braun, Arm=="EVEROLIMUS"),
                         grouping ="set range",
                         ranges = cutpoint$cutpoint$cutpoint,
                         gene = "COL7A1",
                         use_outliers = F,
                         histogram_bins = 5,
                         group_names = factor(c("1:Low","2:High")),
                         histogram_limits_x = c(5,17), 
                         histogram_limits_y = c(0,1), 
                         survival_limit_x=c(0,10))

svg(file = "Suplementary_figure_2/Stage4_KM_curve_braun_NIVOLUMAB.svg")
braun_km_niv[[2]]
dev.off()

svg(file = "Suplementary_figure_2/Stage4_KM_curve_braun_hist_NIVOLUMAB.svg")
braun_km_niv[[1]]
dev.off()

svg(file = "Suplementary_figure_2/Stage4_KM_curve_braun_EVEROLIMUS.svg")
braun_km_evr[[2]]
dev.off()

svg(file = "Suplementary_figure_2/Stage4_KM_curve_braun_hist_EVEROLIMUS.svg")
braun_km_evr[[1]]
dev.off()

