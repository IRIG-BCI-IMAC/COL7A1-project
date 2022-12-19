# suplementary figure 4: correlation between COL7A1 and genes in certain pathway groups

# cox model performed on all expressed genes
library(tidyverse)
library(survival)
library(rstatix)
library(ggpubr)

# Data loading ------------------------------------------------------------
load("0_data/Braun_data_extracted.RData")
load("0_data/EMTAB1980_data_extracted.RData")
load("0_data/GSE167093_data_extracted.RData")
load("0_data/TCGA_data_extracted.RData")

#in case I need symbols for final table
ensembl_to_symbol <- read_tsv(file = "GDCdata/TCGA-KIRC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/000d81dd-9ba4-4852-9090-2bf22f6483f0/09755ce8-ed89-411a-a42f-b3edc4e41eeb.rna_seq.augmented_star_gene_counts.tsv", comment = "N_", skip = 1)[,c(1,2,3)]
vect_names<-ensembl_to_symbol$gene_id
names(vect_names)<-ensembl_to_symbol$gene_name

# PLOTING CORRELATION BETWEEN col7a1 AND GENES INVOLVED IN PROLIFERATION
Prolif_tcga <- expression_tcga %>% 
  select(patient, any_of(c(vect_names[c("PLK1","AURKA","BIRC5","TOP2A","COL7A1")]))) %>% 
  mutate_at(2:6, .funs = ~log2(.+1)) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="TCGA")

Prolif_EMTAB1980 <- EMTAB1980_expression %>% 
  select(patient, any_of(c("PLK1","AURKA","BIRC5","TOP2A","COL7A1"))) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="E-MTAB-1980")

Prolif_GSE167093 <- GSE167093_expression %>% 
  select(patient, any_of(c("PLK1","AURKA","BIRC5","TOP2A","COL7A1"))) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="GSE167093")

Prolif_Braun <- expression_braun %>% 
  select(patient, any_of(c("PLK1","AURKA","BIRC5","TOP2A","COL7A1"))) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="Braun")

proliferation_genes <- rbind(Prolif_tcga, Prolif_EMTAB1980, Prolif_GSE167093, Prolif_Braun)
proliferation_genes <- proliferation_genes %>% 
  mutate(dataset = factor(dataset,
                          levels = c("TCGA",
                                     "E-MTAB-1980",
                                     "GSE167093",
                                     "Braun")))


if(!dir.exists(paths = "Suplementary_figure_4")){
  dir.create(path = "Suplementary_figure_4", recursive = T)
}

svg(file = "Suplementary_figure_4/proliferation_correlation_scaled.svg", height = 15, width = 15)
ggplot(data=proliferation_genes, aes(x=COL7A1 , y=expression))+
  geom_point(size=1)+
  geom_smooth(method = "lm", se = F, linewidth=2, )+
  theme_minimal()+
  stat_cor(label.x = 0, p.accuracy = 0.001, r.accuracy = 0.01, position=position_nudge(y = 1))+
  facet_grid(gene~dataset, scales = "free")+
  geom_rug(col=rgb(.5,0,0,alpha=.2))+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20))
dev.off()


# epithelial mesenchymal transition ----------------------------------------

# PLOTING CORRELATION BETWEEN col7a1 AND GENES INVOLVED IN Metabolism
metabolism_tcga <- expression_tcga %>% 
  select(patient, any_of(c(vect_names[c("ACADL","PDK4","SLC27A2","HAO2","COL7A1")]))) %>% 
  mutate_at(2:6, .funs = ~log2(.+1)) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="TCGA")

metabolism_EMTAB1980 <- EMTAB1980_expression %>% 
  select(patient, any_of(c("ACADL","PDK4","SLC27A2","HAO2","COL7A1"))) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="E-MTAB-1980")

metabolism_GSE167093 <- GSE167093_expression %>% 
  select(patient, any_of(c("ACADL","PDK4","SLC27A2","HAO2","COL7A1"))) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="GSE167093")

metabolism_Braun <- expression_braun %>% 
  select(patient, any_of(c("ACADL","PDK4","SLC27A2","HAO2","COL7A1"))) %>% 
  mutate_at(2:6, .funs = scale) %>% 
  pivot_longer(cols = 2:5, names_to="gene", values_to="expression" ) %>% 
  mutate(dataset="Braun")

metabolism_genes <- rbind(metabolism_tcga, metabolism_EMTAB1980, metabolism_GSE167093, metabolism_Braun)
metabolism_genes <- metabolism_genes %>% 
  mutate(dataset = factor(dataset,
                          levels = c("TCGA",
                                     "E-MTAB-1980",
                                     "GSE167093",
                                     "Braun")))


if(!dir.exists(paths = "Suplementary_figure_4")){
  dir.create(path = "Suplementary_figure_4", recursive = T)
}

svg(file = "Suplementary_figure_4/metabolism_correlation_scaled.svg", height = 15, width = 15)
ggplot(data=metabolism_genes, aes(x=COL7A1 , y=expression))+
  geom_point(size=1)+
  geom_smooth(method = "lm", se = F, linewidth=2, )+
  theme_minimal()+
  stat_cor(label.x = 0, p.accuracy = 0.001, r.accuracy = 0.01, position=position_nudge(y = 1))+
  facet_grid(gene~dataset, scales = "free")+
  geom_rug(col=rgb(.5,0,0,alpha=.2))+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20))
dev.off()

