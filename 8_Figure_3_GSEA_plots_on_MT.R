# Ploting of certain GSEA curves that will be part of figure 3

library(clusterProfiler)
library(enrichplot)

if(!dir.exists(paths = "Figure_3")){
  dir.create(path = "Figure_3", recursive = T)
}


# Gene Position -----------------------------------------------------------
#TCGA data
GP_TCGA <- readRDS("RDS/GSEA_MsigDB_TCGA/Gene_position/c1.all.v2022.1.Hs.symbols.RDS")
svg(file = "Figure_3/Gene_Position_TCGA_MT.svg", width = 10)
gseaplot2(x = GP_TCGA,
          geneSetID = "MT",
          pvalue_table = T, 
          title = "TCGA")
dev.off()


#same for Braun data
GP_BRAUN <- readRDS("RDS/GSEA_MsigDB_BRAUN/Gene_position/c1.all.v2022.1.Hs.symbols.RDS")
svg(file = "Figure_3/Gene_Position_BRAUN_MT.svg", width = 10)
gseaplot2(x = GP_BRAUN,
          geneSetID = "MT",
          pvalue_table = T, 
          title = "BRAUN")
dev.off()

# mitochondrion genes -----------------------------------------------------

#### correlation with of Mitochondria genes
GO_CC_EMTAB1980 <- readRDS("RDS/GSEA_MsigDB_EMTAB1980/GO_CC/c5.go.cc.v2022.1.Hs.symbols.RDS")
svg(file = "Figure_3/Mitochondrion_genes_EMTAB1980.svg", width = 10)
gseaplot2(x = GO_CC_EMTAB1980,
          geneSetID = "GOCC_MITOCHONDRION",
          pvalue_table = T, 
          title = "EMTAB1980")
dev.off()

GO_CC_GSE167093 <- readRDS("RDS/GSEA_MsigDB_GSE167093/GO_CC/c5.go.cc.v2022.1.Hs.symbols.RDS")
svg(file = "Figure_3/Mitochondrion_genes_GSE167093.svg", width = 10)
gseaplot2(x = GO_CC_EMTAB1980,
          geneSetID = "GOCC_MITOCHONDRION",
          pvalue_table = T, 
          title = "GSE167093")
dev.off()

