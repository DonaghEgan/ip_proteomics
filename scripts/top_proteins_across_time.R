top_proteins <- readRDS("/home/degan/ip_proteomics/inputs/protein_matrix_top100.Rds")

gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(top_proteins)),]
top_proteins <- cbind(top_proteins, gene_name_matrix[, "Gene.names", drop=FALSE])
top_proteins <- top_proteins[!duplicated(top_proteins[,"Gene.names"]),]
rownames(top_proteins) <- NULL
top_proteins <- column_to_rownames(top_proteins, "Gene.names")

cor_proteins <- cor(t(top_proteins), method = "spearman")

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/PDCD1_module_heatmap.pdf", height = 5, width = 5)
protein_heatmap <- ComplexHeatmap::Heatmap(cor_proteins,name = "Corr", show_row_names = T,
                                           show_column_names = FALSE, 
                                           row_names_gp = gpar(fontsize = 4))
dev.off()

## gsea on gene modules from correlation plot ####
cluster_1 <- c("PDCD1", "PTPN11", "HSPA1A;HSPA1B", "USF1", "LGALS9", "MSN")

cluster_2 <- c("RAC1;RAC2","HSP90AB1", "GAPDH", "PFN1", "CCT4","ATP6V1E1","EEF1D",
               "HSP90B1", "ATP5A1", "PKM", "HMGB1;HMGB1P1", "CCT3", "EEF1A1;EEF1A1P5")

cluster_5 <- c("TCF7","SLC16A1", "TBL2", "SMN1", "GRB2", "CDKN2AIP", "ARHGEF1",
               "ATR", "MYO18A", "ATRIP")



## fgsea ####
ensembl.pathway <- sbgn.gsets(id.type = "SYMBOL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = T, #T
                              #database = "MetaCyc", 
                              truncate.name.length = 100)

msigdbr_df <- msigdbr(species = "human", category = "C3")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

genes <- diff_exp_protein %>% pull(logFC, Gene.names)
genes <- genes[cluster_5]

fgseaRes <- fgsea(pathways = ensembl.pathway, stats = genes)
fgseaRes$pathway <- gsub(".*::","",fgseaRes$pathway)
fgseaRes <- fgseaRes %>% arrange(desc(abs(NES)))

## plotting expression of 5 genes common to all time points across time ####
common_genes <- diff_exp_protein[which(diff_exp_protein$logFC >=1.5 & diff_exp_protein$adj.P.Val <0.05),]
common_genes <- data.frame(table(common_genes$Gene.names))
common_genes <- common_genes[which(common_genes$Freq >= 4), ]
names(common_genes)[1] <- "Gene.names"

diff_exp_protein_common <- diff_exp_protein[which(diff_exp_protein$Gene.names %in% common_genes$Gene.names),]
common_genes <- cbind(diff_exp_protein_common[, c("rn","Gene.names")], common_genes)
common_genes <- common_genes[!duplicated(common_genes)]


## imputed protein matrix 
imputed_proteins <- readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")

imputed_common_genes <- imputed_proteins[which(rownames(imputed_proteins) %in% common_genes$rn),]
imputed_common_genes_pd1 <- data.frame(imputed_common_genes[,which(colnames(imputed_common_genes) %in% rownames(pd1_only))])
imputed_common_genes_pd1$rn <- rownames(imputed_common_genes_pd1)
imputed_common_genes_pd1 <- data.frame(cbind(imputed_common_genes_pd1, common_genes[,c(1,2)]))
rownames(imputed_common_genes_pd1) <- NULL
imputed_common_genes_pd1 <- column_to_rownames(imputed_common_genes_pd1, "Gene.names")
imputed_common_genes_pd1 <- imputed_common_genes_pd1[,c(1:96)]

imputed_common_genes_pd1 <- data.frame(t(imputed_common_genes_pd1))
imputed_common_genes_pd1$time <- sapply(strsplit(rownames(imputed_common_genes_pd1), "_"), "[", 2)

imputed_common_genes_pd1$time <- gsub("min", "", imputed_common_genes_pd1$time) 
imputed_common_genes_pd1$time <- gsub("hr", "", imputed_common_genes_pd1$time) 
imputed_common_genes_pd1$time <- as.numeric(imputed_common_genes_pd1$time)
imputed_common_genes_pd1 <- aggregate(. ~ time, imputed_common_genes_pd1, mean)

imputed_common_genes_pd1 <- melt(imputed_common_genes_pd1, id.vars = "time")

pdf("/home/degan/ip_proteomics/figures/common_genes_overtime.pdf", width = 4.5, height = 3)
ggplot(data=imputed_common_genes_pd1, aes(x=time, y=value, group=variable)) +
  geom_line(aes(color=variable), linetype = "dashed")+
  geom_point(aes(color=variable)) + ylab("Mean Abundance") + theme_bw() +
  scale_color_brewer(palette="Set2")
dev.off()







