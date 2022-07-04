################################################################################
# Donagh Egan
# POIAZ 
# Date: 30nd June 2022
# 

# Description: Analyze clusters from ip proteomics in the whole and phospho data 
################################################################################
library(ggsci)

## load in data ####
################################################################################

cluster_assignment <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
w_proteomics <- readRDS("/mnt/Data/Vadim/POIAZ/Donagh/proteomics22_unified/processing_&_normalisation/saved_objects/processed_proteomics.Rds")
exp_w_proteomics <- readRDS("/mnt/Data/Vadim/POIAZ/Donagh/proteomics22_unified/processing_&_normalisation/saved_objects/experimental_design.Rds")
gene_meta_data <-  readRDS("/mnt/Data/Vadim/POIAZ/Donagh/proteomics22_unified/processing_&_normalisation/saved_objects/phospho_gene_names.Rds")

## dividing exp design ctrl and pd1 ####
################################################################################

exp_pd1 <- exp_w_proteomics[exp_w_proteomics$type == "PD1",]
exp_ctrl <- exp_w_proteomics[exp_w_proteomics$type == "CNTR",]

## genes in PD1 ####
################################################################################

cluster_assignment <- cluster_assignment[cluster_assignment$cluster == 13,]
pd1_genes_meta <- gene_meta_data[gene_meta_data$Gene.names %in% cluster_assignment$Gene.names,]

## Dividing counts into pd1 and ctrl ####
pd1_proteins <- cluster_assignment[cluster_assignment$cluster == 13,]
diff_exp_protein <- diff_exp_protein[diff_exp_protein$type == "PD1vsCtr",]
pd1_up <- diff_exp_protein[diff_exp_protein$logFC > .5 & diff_exp_protein$adj.P.Val < 0.05, ]
pd_up_cluster <- intersect(pd1_proteins$Row.names, pd1_up$rn)

pd_up_cluster <- strsplit(pd_up_cluster, ";", fixed=F)
pd_up_cluster <- Reduce(c,pd_up_cluster)

pd1_proteomics <- w_proteomics[rownames(w_proteomics) %in% pd_up_cluster,]
pd1_ctrl <- pd1_proteomics[,colnames(pd1_proteomics) %in% rownames(exp_ctrl)]
pd1_only <- pd1_proteomics[,colnames(pd1_proteomics) %in% rownames(exp_pd1)]



ctrl_pd1_avg <- melt(data.frame(pd1=rowMeans(pd1_only), ctrl=ctrl_avg <- rowMeans(pd1_ctrl)))

ggplot(ctrl_pd1_avg, aes(x = value, fill = variable)) + geom_density(alpha = 0.5) +
  theme_bw() + xlim(10,30) + ylim(0, 0.25) + scale_fill_npg()

## testing for phospho ####
################################################################################
phospho_proteomics <- readRDS("/mnt/Data/Vadim/POIAZ/Donagh/proteomics22_unified/processing_&_normalisation/saved_objects/imputed_normalised_phospho.Rds")
phospho_proteomics_assay <- data.frame(assay(phospho_proteomics))
colnames(phospho_proteomics_assay) <- as.numeric(gsub("X", "", colnames(phospho_proteomics_assay)))
phospho_proteomics_assay$gene_name <- phospho_proteomics@GeneSymbol
phospho_proteomics_assay <- phospho_proteomics_assay[phospho_proteomics_assay$gene_name %in% pd1_genes_meta$Gene.names,]


conditions_ctrl <- conditions[conditions$type == "CNTR",]
conditions_pd1 <- conditions[conditions$type == "PD1",]

phospho_proteomics_ctrl <- phospho_proteomics_assay[, colnames(phospho_proteomics_assay) %in% rownames(conditions_ctrl)]
phospho_proteomics_pd1 <- phospho_proteomics_assay[, colnames(phospho_proteomics_assay) %in% rownames(conditions_pd1)]

ctrl_pd1_avg_phos <- melt(data.frame(pd1=rowMedians(as.matrix(phospho_proteomics_pd1)), ctrl=rowMedians(as.matrix(phospho_proteomics_ctrl)))

ggplot(ctrl_pd1_avg_phos, aes(x = value, fill = variable)) + geom_density(alpha = 0.5) +
  theme_bw() + scale_fill_npg() + xlim(-5, 5)








