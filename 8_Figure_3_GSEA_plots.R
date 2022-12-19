# Making a barplot of GSEA results
library(tidyverse)

halmark_EMTAB1980 <- read_csv("CSV/GSEA_MsigDB_EMTAB1980/Halmark/h.all.v2022.1.Hs.symbols.csv") %>% 
  mutate(dataset='E-MTAB-1980') 
halmark_GSE167093 <- read_csv("CSV/GSEA_MsigDB_GSE167093/Halmark/h.all.v2022.1.Hs.symbols.csv") %>% 
  mutate(dataset='GSE167093')
halmark_TCGA <- read_csv("CSV/GSEA_MsigDB_TCGA/Halmark/h.all.v2022.1.Hs.symbols.csv") %>% 
  mutate(dataset='TCGA')
halmark_Braun <- read_csv("CSV/GSEA_MsigDB_BRAUN/Halmark/h.all.v2022.1.Hs.symbols.csv") %>% 
  mutate(dataset='Braun et al.')

TCGA_select <- head(halmark_TCGA, n = 10)
EMTAB1980_select <- head(halmark_EMTAB1980, n=10)
GSE167093_select <- head(halmark_GSE167093, n=10)
Braun_select <- head(halmark_Braun, n=10)

# selection vector
selection_vector <- c(TCGA_select$ID,EMTAB1980_select$ID,GSE167093_select$ID, Braun_select$ID )
selection_vector <- unique(selection_vector)

joined <- rbind(halmark_EMTAB1980,halmark_GSE167093,halmark_TCGA,halmark_Braun) %>% 
  filter(ID %in% selection_vector) %>% 
  mutate(dataset=factor(dataset, levels = c("TCGA","E-MTAB-1980","GSE167093","Braun et al.")))%>% 
  mutate(significant= ifelse(p.adjust<0.05, "FDR < 0.05", "FDR > 0.05")) %>% 
  mutate(enrichment= ifelse(NES>0, "Positive", "Negative")) %>% 
  mutate(Significance= paste0(significant, ", ", enrichment)) %>% 
  mutate(ID=str_replace_all(str_split_fixed(ID, "_",2)[,2], "_", " "))


# making a barplot 

if(!dir.exists(paths = "Figure_3")){
  dir.create(path = "Figure_3", recursive = T)
}

svg(file = 'Figure_3/HALMARK_pathways_fig_1.svg', width =12, height = 8)
ggplot(data = joined, mapping = aes(x = NES, y=reorder(ID,NES),
                                    fill=Significance))+
  scale_fill_manual(values = c("FDR < 0.05, Positive"="red",
                               "FDR > 0.05, Positive"="orange",
                               "FDR < 0.05, Negative"="blue",
                               "FDR > 0.05, Negative"="skyblue"))+
  scale_color_manual(values = c("FDR < 0.05, Positive"="red",
                                "FDR > 0.05, Positive"="orange",
                                "FDR < 0.05, Negative"="blue",
                                "FDR > 0.05, Negative"="skyblue"))+
  geom_col(position = "dodge", width = 0.1, show.legend = T,
           color="BLACK")+
  geom_point(aes(color=Significance), size=3 ,show.legend = F)+
  geom_text(mapping = aes(label=round(NES, 1), x = (-0.8 * sign(NES))))+
  geom_vline(xintercept = 0)+
  theme_minimal()+
  facet_grid(~dataset)
dev.off()

