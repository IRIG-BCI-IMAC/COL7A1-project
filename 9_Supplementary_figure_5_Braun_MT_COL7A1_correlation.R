# cox model performed on all expressed genes
library(tidyverse)
library(survival)
library(rstatix)
library(ggpubr)

if(!dir.exists(paths = "Suplementary_figure_5")){
  dir.create(path = "Suplementary_figure_5", recursive = T)
}

# Data loading ------------------------------------------------------------
load("0_data/Braun_data_extracted.RData")

# Data transformation -----------------------------------------------------

# I will compare MT genes expression acros stage, grade, and COL7A1 levels

mt_genes <- colnames(expression_braun) 
mt_genes <- c(mt_genes[str_sub(mt_genes,1, 3)=="MT-"], "COL7A1")

mitochondrial_data <- expression_braun %>% 
  select(patient, any_of(c(mt_genes, "COL7A1"))) %>% 
  mutate_at(2:15, .funs = scale)

mitochondrial_data <- mitochondrial_data %>% 
  pivot_longer(cols = 2:14, names_to="gene", values_to="expression" )

svg(file = "Suplementary_figure_5/correlation_with_MT_genes_Braun.svg", height = 15, width = 15)
ggplot(data=mitochondrial_data, aes(x=COL7A1 , y=expression))+
  geom_point(size=1)+
  geom_smooth(method = "lm", se = F, linewidth=2, )+
  theme_minimal()+
  stat_cor(label.x = 0, p.accuracy = 0.001, r.accuracy = 0.01, position=position_nudge(y = 1))+
  facet_wrap(~gene, scales = "free")+
  geom_rug(col=rgb(.5,0,0,alpha=.2))+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 20))

dev.off()
