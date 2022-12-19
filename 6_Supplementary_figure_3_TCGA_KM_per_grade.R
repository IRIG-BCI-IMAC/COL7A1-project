# cox model performed on all expressed "COL7A1"s
library(tidyverse)
library(survival)

load("0_data/TCGA_data_extracted.RData")
source(file = "1_functions/kaplan_meier_curve.R")

if(!dir.exists(paths = "Suplementary_figure_3/")){
  dir.create(path = "Suplementary_figure_3/", recursive = T)
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

# since G1 has too few patients, I will combine G1 and G2
clinical_tcga<-clinical_tcga %>% 
  mutate(Histologic_grade=as.character(Histologic_grade)) %>% 
  mutate(Histologic_grade=ifelse(Histologic_grade %in% c("G1","G2"),"G1/2",Histologic_grade)) %>% 
  mutate(Histologic_grade = factor(Histologic_grade)) %>% 
  filter(!Histologic_grade %in% c("","GX"))

## Histologic_grade I+II
Histologic_grade_I_II <- clinical_tcga %>% 
  filter(Histologic_grade=="G1/2")
Histologic_grade_I_II<-expression_tcga %>% 
  select(patient, "COL7A1") %>% 
  semi_join(Histologic_grade_I_II)


KM_TCGA_COL7A1_Histologic_grade_1_2 <- KM_curve(mRNA_counts = Histologic_grade_I_II, clinical_data = clinical_tcga, gene = "COL7A1",
                                                grouping = "set range", censoring = T, use_outliers = T,
                                                group_names=c("Low","High"),
                                                ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,5), 
                                                histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))
## Histologic_grade III

Histologic_grade_III <- clinical_tcga %>% 
  filter(Histologic_grade=="G3")
Histologic_grade_III<-expression_tcga %>% 
  select(patient, "COL7A1") %>% 
  semi_join(Histologic_grade_III)

KM_TCGA_COL7A1_Histologic_grade_3 <- KM_curve(mRNA_counts = Histologic_grade_III, clinical_data = clinical_tcga, gene = "COL7A1",
                                   grouping = "set range", censoring = T, use_outliers = T,
                                   group_names=c("Low","High"),
                                   ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,5), 
                                   histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))

## Histologic_grade IV

Histologic_grade_IV <- clinical_tcga %>% 
  filter(Histologic_grade=="G4")
Histologic_grade_IV<-expression_tcga %>% 
  select(patient, "COL7A1") %>% 
  semi_join(Histologic_grade_IV)

KM_TCGA_COL7A1_Histologic_grade_4 <- KM_curve(mRNA_counts = Histologic_grade_IV, clinical_data = clinical_tcga, gene = "COL7A1",
                                   grouping = "set range", censoring = T, use_outliers = T,
                                   group_names=c("Low","High"),
                                   ranges=cutpoint$cutpoint$cutpoint,histogram_limits_x = c(0,5), 
                                   histogram_limits_y = c(0,50),  survival_limit_x=c(0,10))

# histogram
hist_data<-clinical_tcga %>% filter(Histologic_grade!="")
hist_data <- expression_tcga %>% 
  semi_join(hist_data) 

histogram <- KM_curve(mRNA_counts = hist_data, clinical_data = clinical_tcga, gene = "COL7A1",
                      grouping = "set range", censoring = T, use_outliers = T,
                      group_names=c("Low","High"), histogram_bins=2,
                      ranges=cutpoint$cutpoint$cutpoint,facet_by = "Histologic_grade",
                      histogram_limits_x = c(0,7.5), 
                      histogram_limits_y = c(0,1),  survival_limit_x=c(0,10))


svg("Suplementary_figure_3/Histologic_grade__1_2.svg")
KM_TCGA_COL7A1_Histologic_grade_1_2[[2]]+ggtitle("G1/G2")
dev.off()

svg("Suplementary_figure_3/Histologic_grade__3.svg")
KM_TCGA_COL7A1_Histologic_grade_3[[2]]+ggtitle("G3")
dev.off()

svg("Suplementary_figure_3/Histologic_grade__4.svg")
KM_TCGA_COL7A1_Histologic_grade_4[[2]]+ggtitle("G4")
dev.off()

svg("Suplementary_figure_3/Histologic_grade__histogram.svg", width = 24)
histogram[[1]]
dev.off()



