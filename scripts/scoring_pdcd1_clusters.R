library(GSVA)

gene_name_matrix <- pgroups.meta[which(pgroups.meta$Protein.IDs %in% rownames(imputed_nostim)),]
imputed_nostim_genenames <- merge(imputed_nostim, gene_name_matrix, by.x = 0, by.y = "Protein.IDs")
imputed_nostim_genenames <-  imputed_nostim_genenames[!duplicated(imputed_nostim_genenames[,"Gene.names"]),]
rownames(imputed_nostim_genenames) <- NULL
imputed_nostim_genenames <- imputed_nostim_genenames %>% column_to_rownames("Gene.names")
imputed_nostim_genenames$Row.names <- NULL

cluster_1 <- c("PDCD1", "PTPN11", "HSPA1A;HSPA1B", "USF1", "LGALS9", "MSN")

cluster_2 <- c("RAC1;RAC2","HSP90AB1", "GAPDH", "PFN1", "CCT4","ATP6V1E1","EEF1D",
               "HSP90B1", "ATP5A1", "PKM", "HMGB1;HMGB1P1", "CCT3", "EEF1A1;EEF1A1P5")

cluster_4 <- c("SLC25A3", "SLC25A6", "HSPA6", "HSPA8", "ZRANB2", "C3", "TYMS", 
               "YWHAG", "YWHAQ")

cluster_3 <- c("MTF2", "HMGN2", "LCK", "HNRNPU", "HIST1H1B", "HIST1H1E", "HIST1H1C",
               "CAPZA1", "FAM133B;FAM133A", "HIST1H3A;HIST3H3;H3F3C", "MYH14")

cluster_5 <- c("ZAP70", "CD3E", "CD247")

cluster_6 <- c("TCF7","SLC16A1", "TBL2", "SMN1", "GRB2", "CDKN2AIP", "ARHGEF1",
               "ATR", "MYO18A", "ATRIP")


pd1_modules <- list("cluster 1" = cluster_1, "cluster 2" = cluster_2, 
                    "cluster 3" = cluster_3, "cluster 4" = cluster_4, 
                    "cluster 5" = cluster_5, "cluster 6" = cluster_6)

pd1_module_scores <- gsva(as.matrix(imputed_nostim_genenames), pd1_modules, mx.diff=FALSE, 
                          verbose=TRUE, parallel.sz=1, method = "ssgsea")

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/pd1_modules_scored.pdf", height = 4, width = 5)
pheatmap::pheatmap(pd1_module_scores, scale = "row", show_colnames = FALSE,
                   annotation = annotation_nostim[,c(1,2,3)], treeheight_row = 15,
                   treeheight_col = 25)
dev.off()

## Note: clustering performed on rows as well. Cluster 2 and 3 are specific to
## pembrolizumab and Csa. Cluster 5 more specific to nivolumab. Cluster 1 and
## 4 are antibody agnostic. 
# cluster 6 early timepoints 


