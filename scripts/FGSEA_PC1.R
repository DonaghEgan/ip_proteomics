library(readxl)
library(SBGNview)
library(fgsea)
library(msigdbr)
library(openxlsx)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(GO.db)
library(xtable)
library(annotate)
library(genefilter)
library(GOstats)

## Setting up pathways ####
msigdbr_df <- msigdbr(species = "human", category = "C2", subcategory = "CP")
pathways = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

## Run PCA script prior to this script to load in variables ####
rownames(loadings_top_pc1) <- sub(";.*", "", rownames(loadings_top_pc1))

loadings_pc1 <- PCAloadings %>% select(PC1)

gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(loadings_pc1)),]
loadings_pc1 <- merge(gene_name_matrix[,c("Protein.IDs", "Gene.names")], loadings_pc1, by.x = "Protein.IDs", by.y = 0)

genes <- loadings_pc1 %>% pull(PC1, Gene.names)

fgseaRes <- fgsea(pathways = pathways, stats = genes, nPermSimple = 10000)










haveGo <- sapply(mget(genes_ID$ENTREZID, org.Hs.egGO),
                  function(x) {
                    if (length(x) == 1 && is.na(x))
                      FALSE
                    else TRUE
                    })

numNoGO <- sum(!haveGo)
genes_ID <- genes_ID[haveGo, ]

## Define gene universe based on results of non-specific filtering ####

universe <- genes_ID$ENTREZID
top_pc1 <-  PCAloadings %>% arrange(desc(abs(PC1))) %>% select(PC1)
rownames(top_pc1) <- sub(";.*", "", rownames(top_pc1))
pc1_go <- GO.df[which(GO.df$UNIPROTKB %in% rownames(top_pc1)),]
pc1_go <- pc1_go[!duplicated(pc1_go$UNIPROTKB), ]

test <- merge(pc1_go, top_pc1, by.x="UNIPROTKB", by.y=0)
test <- merge(GO.annotation, test, by.x="GOID",by.y="GOID")

## order samples based on PC1 value ###
pc1_order <- PCAvalues %>% arrange(PC1) %>% select(PC1)

## protein exp ####
pc1_data <- imputed_data[which(rownames(imputed_data) %in% test$UNIPROTKB),]
pc1_data <- pc1_data[rownames(pc1_order)]
pc1_data <- data.frame(t(pc1_data))
  
pc1_data <- merge(pc1_data, test[, c(2,4)], by.x=0,by.y="UNIPROTKB")
pc1_data <- column_to_rownames(pc1_data, "Row.names")
pc1_data <- data.frame(t(pc1_data))
pc1_data <- scale(pc1_data)


col_ha <- ComplexHeatmap::HeatmapAnnotation("GOID" = data.frame(test[,4]), which = "col")

row_ha <- ComplexHeatmap::rowAnnotation(sample_id = rownames(pc1_data)[1:nrow(pc1_data)-1])

ComplexHeatmap::Heatmap(as.matrix(pc1_data), cluster_rows = F)








