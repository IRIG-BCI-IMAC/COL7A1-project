# cox model performed on all expressed genes
library(tidyverse)
library(survival)
library(rstatix)
library(ggpubr)

if(!dir.exists(paths = "Suplementary_figure_5")){
  dir.create(path = "Suplementary_figure_5", recursive = T)
}

# Data loading ------------------------------------------------------------
load("0_data/TCGA_data_extracted.RData")

#in case I need symbols for final table
ensembl_to_symbol <- read_tsv(file = "GDCdata/TCGA-KIRC/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/000d81dd-9ba4-4852-9090-2bf22f6483f0/09755ce8-ed89-411a-a42f-b3edc4e41eeb.rna_seq.augmented_star_gene_counts.tsv", comment = "N_", skip = 1)[,c(1,2,3)]
vect_names<-ensembl_to_symbol$gene_id
names(vect_names)<-ensembl_to_symbol$gene_name

# Data transformation -----------------------------------------------------

# I will compare MT genes expression acros stage, grade, and COL7A1 levels

mt_genes <- ensembl_to_symbol %>% 
  filter(str_sub(gene_name,1, 3)=="MT-")
vect_names <- vect_names[vect_names %in% c(mt_genes$gene_id, "ENSG00000114270.17")]

mitochondrial_data <- expression_tcga %>% 
  select(patient, any_of(c(mt_genes$gene_id, "ENSG00000114270.17"))) %>% 
  mutate_at(2:15, .funs = ~log2(.+1)) %>% 
  mutate_at(2:15, .funs = scale) 

vect_names <- vect_names[vect_names %in% colnames(mitochondrial_data)]

mitochondrial_data <- mitochondrial_data%>% 
  rename(vect_names) %>% 
  pivot_longer(cols = 2:14, names_to="gene", values_to="expression" )

svg(file = "Suplementary_figure_5/correlation_with_MT_genes_TCGA.svg", height = 15, width = 15)
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
