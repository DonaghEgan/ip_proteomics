################################################################################
# Donagh Egan
# POIAZ 
# Date: 16th March 2022
# 

# Description: Read in all (non-filtered), differentially bound proteins from limma model
# and run FGSEA on for each comparison made.
# Uses output from 
################################################################################

## Library #### 
################################################################################

library(readxl)
library(SBGNview)
library(fgsea)
library(msigdbr)
library(openxlsx)
library(topGO)
library(clusterProfiler)

## Read in differentially bound proteins ####
diff_bound <- read_excel("/home/degan/ip_proteomics/inputs/differentially_bound_proteins_long.xlsx")
processed_ip <- readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")


## Setting up pathways ####
msigdbr_df <- msigdbr(species = "human", category = "C7")
pathways = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

## SBGN view for each comparison ####
data("pathways.info")
ensembl.pathway <- sbgn.gsets(id.type = "SYMBOL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = T, 
                              #database = "MetaCyc", 
                              truncate.name.length = 100)

## FGSEA for each comparison ####
fgsea_result <- list()
for (condition in unique(diff_bound$type)) {
     proteins.condition <- diff_bound[diff_bound$type == condition,] %>% arrange(desc(abs(t)))
     genes <- proteins.condition %>% pull(t, Gene.names)
     names(genes) <- sub(";.*", "", names(genes))
     fgseaRes <- fgsea(pathways = pathways, stats = genes, nPermSimple = 10000)
     fgseaRes <- fgseaRes[which(fgseaRes$padj <0.1),]
     
     fgsea_result[[condition]] = fgseaRes
}

sbgn_result <- list()
for (condition in unique(diff_bound$type)) {
     proteins.condition <- diff_bound[diff_bound$type == condition,] %>% arrange(desc(abs(t)))
     genes <- proteins.condition %>% pull(t, Gene.names)
     names(genes) <- sub(";.*", "", names(genes))
     fgseaRes <- fgsea(pathways = ensembl.pathway, stats = genes, nPermSimple = 10000)
     fgseaRes <- fgseaRes[which(fgseaRes$padj <0.05),]
     sbgn_result[[condition]] = fgseaRes
}

## GO classification ####
geneList <- diff_bound$adj.P.Val
names(geneList) <- diff_bound$Gene.names
names(geneList) <- sub(";.*", "", names(geneList))
names(geneList) <- na.omit(names(geneList))

gene <- names(geneList)[abs(geneList) > 2]


gene <- AnnotationDbi::select(org.Hs.eg.db,keys = gene, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")

go_result <- list()
for (condition in unique(diff_bound$type)) {
  proteins.condition <- diff_bound[diff_bound$type == condition, c("adj.P.Val", "logFC", "Gene.names")] 
  genes <- proteins.condition[which(proteins.condition$logFC > 0 & proteins.condition$adj.P.Val <0.1), ]
  genes <- genes$Gene.names
  genes <- AnnotationDbi::select(org.Hs.eg.db,keys = genes, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
  ggo <- enrichGO(gene = genes$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE, pvalueCutoff  = 0.01, qvalueCutoff  = 0.01, )
  go_result[[condition]] = ggo
  
}

## saving results ####
write.xlsx(fgsea_result, file = "/home/degan/ip_proteomics/inputs/pathway_enrichment/H_fgsea_res_ip_long.xlsx")
write.xlsx(fgsea_result, file = "/home/degan/ip_proteomics/inputs/pathway_enrichment/C7_fgsea_res_ip_long.xlsx")
write.xlsx(sbgn_result, file = "/home/degan/ip_proteomics/inputs/pathway_enrichment/SBGN_fgsea_res_ip_long.xlsx")
write.xlsx(go_result, file = "/home/degan/ip_proteomics/inputs/pathway_enrichment/Go_enrichment_res_ip_long.xlsx")











