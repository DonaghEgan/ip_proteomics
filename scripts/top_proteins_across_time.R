top_proteins <- readRDS("/home/degan/ip_proteomics/inputs/protein_matrix_top100.Rds")

gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(top_proteins)),]
top_proteins <- cbind(top_proteins, gene_name_matrix[, "Gene.names", drop=FALSE])
top_proteins <- top_proteins[!duplicated(top_proteins[,"Gene.names"]),]
rownames(top_proteins) <- NULL
top_proteins <- column_to_rownames(top_proteins, "Gene.names")

## top proteins in PD1 samples only ####
pd1_only <- annotation_nostim[which(substr(annotation_nostim$antibody, 1, 3) == "PD1"),]
top_proteins_pd1 <- top_proteins[,which(colnames(top_proteins) %in% rownames(pd1_only))]

cor_proteins <- cor(t(top_proteins_pd1), method = "spearman")

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/PDCD1_module_heatmap.pdf", height = 5, width = 5)
protein_heatmap <- ComplexHeatmap::Heatmap(cor_proteins,name = "Corr", show_row_names = T,
                                           show_column_names = FALSE, 
                                           row_names_gp = gpar(fontsize = 4))
dev.off()

## gsea on gene modules from correlation plot ####
cluster_1 <- c("TRMT61A","SPOP","CPSF7","LACTB","TCF7","TNRC6C","CACTIN",
               "TRAF3IP3","ALYREF","POLR2L","TIGD2","HSPH1","DNAJA3","MOV10",
               "RBMS1","DDX6","EIF4ENIF1","AGO3","MYO7B","RSBN1","RBMX2",
               "CABIN1","RBBP6")

cluster_2 <- c("SLC25A3","PPIE", "ABCF2", "VPS35", "HNRNPK", "ERCC3", "SGPL1", "ZNF830",
               "KPNA1", "UBN1", "RSRC1", "SREK1IP1", "RAC1", "RAC2","GPATCH4", "AHNAK",
               "RNPC3", "TFAP4", "MSANTD4", "LONP1", "SMU1", "DDX42", "WHSC1L1","TIMM13",
               "PTPN11", "CEBPB", "EXOSC5", "IWS1", "SNRNP27", "SFSWAP","CTNNBL1", "SNIP1",
               "ATAD2B", "SCAF1", "ISG20L2", "PDCD1", "C3orf17", "SMARCA2", "CENPV", "SP1",
               "DR1", "MAX", "TIAL1", "USF1", "JPH1", "AP2B1", "SUN1", "URB2", "RBX1", "TAF1D",
               "DIAPH1", "NACA", "SSSCA1", "CASP8AP2", "ATP5C1", "OR10AG1", "AMOTL2", "RREB1",
               "RRP1B", "MDC1")

## fgsea ####
ensembl.pathway <- sbgn.gsets(id.type = "SYMBOL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = T, #T
                              #database = "MetaCyc", 
                              truncate.name.length = 100)

msigdbr_df <- msigdbr(species = "human", category = "C2")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

genes <- diff_exp_protein %>% pull(logFC, Gene.names)
genes <- genes[cluster_1]
genes <- na.omit(genes)

fgseaRes <- fgsea(pathways = pathwaysH, stats = genes)
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







