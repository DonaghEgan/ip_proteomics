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
msigdbr_df <- msigdbr(species = "human", category = "H")
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
     fgseaRes$condition <- condition
     
     fgsea_result[[condition]] = fgseaRes
}

sbgn_result <- list()
for (condition in unique(diff_bound$type)) {
     proteins.condition <- diff_bound[diff_bound$type == condition,] %>% arrange(desc(abs(t)))
     genes <- proteins.condition %>% pull(t, Gene.names)
     names(genes) <- sub(";.*", "", names(genes))
     fgseaRes <- fgsea(pathways = ensembl.pathway, stats = genes, nPermSimple = 10000)
     fgseaRes <- fgseaRes[which(fgseaRes$padj <0.05),]
     fgseaRes$condition <- condition
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
write.xlsx(fgsea_result, file = "/home/degan/ip_proteomics/inputs/pathway_enrichment/Canonical_KEGG_Reactome_res_ip_long.xlsx")
write.xlsx(fgsea_result, file = "/home/degan/ip_proteomics/inputs/pathway_enrichment/C2_res_ip_long.xlsx")


five_t <- fgsea_result[["ST5vNsT0"]]
t24_t <- fgsea_result[["ST24vST20"]]
common_pathways <- intersect(five_t$pathway, t24_t$pathway)
five_t <- five_t[which(five_t$pathway %in% common_pathways),]
t24_t <- t24_t[which(t24_t$pathway %in% common_pathways),]
cytosolic_ribsome <- five_t$leadingEdge[[6]]

## Plotting number of enriched pathways ####
no_pathways.df <- data.frame(t(data.frame(ST5vNsT0 =  nrow(fgsea_result[["ST5vNsT0"]]),
                             ST20vST5 = nrow(fgsea_result[["ST20vST5"]]),
                             ST24vST20 = nrow(fgsea_result[["ST24vST20"]]))))
no_pathways.df$comparison <- rownames(no_pathways.df)
colnames(no_pathways.df)[1] <- "val"

pdf("/home/degan/ip_proteomics/figures/no_enrichedpath_pd-time.pdf", width = 4, height = 4)
positions <- no_pathways.df$comparison
ggplot(no_pathways.df, aes(x=comparison, y=val, fill=comparison)) + ylab("Number of Pathways") +
  geom_bar(stat="identity") + scale_x_discrete(limits = positions) + theme_bw()
dev.off()


## Analysing Hallmarks PD1 vs Ctrl ####
comparisons_hallmarks <- do.call(rbind, fgsea_result)
comparisons_hallmarks <- comparisons_hallmarks[1:19,c("pathway", "NES", "condition")]
comparisons_hallmarks <- reshape(comparisons_hallmarks, idvar = "pathway", timevar = "condition", direction = "wide")
rownames(comparisons_hallmarks) <- NULL
comparisons_hallmarks <- column_to_rownames(comparisons_hallmarks, "pathway")
colnames(comparisons_hallmarks) <- gsub("NES.", "", colnames(comparisons_hallmarks))
comparisons_hallmarks <- comparisons_hallmarks[,-1]

pdf("/home/degan/ip_proteomics/figures/Hallmark_heatmap_ctrvspd1.pdf", height = 4, width = 6)
pheatmap(as.matrix(comparisons_hallmarks), cluster_rows = F, cluster_cols = F)
dev.off()

comparisons_hallmarks_1 <- comparisons_hallmarks[,1, drop=FALSE]
pheatmap(as.matrix(comparisons_hallmarks_1), cluster_rows = F, cluster_cols = F)

## plotting SBGNview results ####
comparisons_SBGN <- do.call(rbind, sbgn_result)
unique_ST24_vs_ST20 <- setdiff(sbgn_result[["ST24vST20"]][["pathway"]], c(sbgn_result[["PD1_24vsCon_24"]][["pathway"]],
                                                                          sbgn_result[["PD1_20vsCon_20"]][["pathway"]],
                                                                          sbgn_result[["PD1_24vsCon_24"]][["pathway"]],
                                                                          sbgn_result[["PD1vsCtr"]][["pathway"]]))

unique_ST20_vs_ST5 <- setdiff(sbgn_result[["ST20vST5"]][["pathway"]], c(sbgn_result[["PD1_24vsCon_24"]][["pathway"]],
                                                                         sbgn_result[["PD1_20vsCon_20"]][["pathway"]],
                                                                         sbgn_result[["PD1_24vsCon_24"]][["pathway"]],
                                                                         sbgn_result[["PD1vsCtr"]][["pathway"]]))

unique_ST5_vs_ST0 <- setdiff(sbgn_result[["ST5vNsT0"]][["pathway"]], c(sbgn_result[["PD1_24vsCon_24"]][["pathway"]],
                                                                        sbgn_result[["PD1_20vsCon_20"]][["pathway"]],
                                                                        sbgn_result[["PD1_5vsCon_5"]][["pathway"]],
                                                                        sbgn_result[["PD1vsCtr"]][["pathway"]]))




## Proteins unique to time ####
unique_pd1_time <- setdiff(diff_bound$Gene.names[diff_bound$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")],
                           diff_bound$Gene.names[diff_bound$type %in% c("PD1_24vsCon_24", "PD1_20vsCon_20", "PD1_5vsCon_5",
                                                                        "PD1vsCtr")])


