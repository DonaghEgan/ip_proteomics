## Library ####
library(limma)
library(data.table)
library(dplyr)
library(stringr)
library(fgsea)
library(msigdbr)
library(SBGNview)
library(SummarizedExperiment)
library(UpSetR)
library(ComplexHeatmap)
library(reshape2)
library(RColorBrewer)
library(RVenn)


## Limma model factors ####
condition <- substr(sapply(strsplit(colnames(imputed_nostim), "_"), "[", 1), 1, 3)
antibody <- substr(sapply(strsplit(colnames(imputed_nostim), "_"), "[", 1), 4,6)
time <-  sapply(strsplit(colnames(imputed_nostim), "_"), "[", 2)

time <- gsub("min", "", time)
time <- gsub("hr", "", time)

condition_time <- paste(condition,time, sep="_")
designlimmaP <- as.data.frame(model.matrix(~ 0 + condition_time + antibody))

colnames(designlimmaP) <- gsub("condition_time", "", colnames(designlimmaP))
colnames(designlimmaP) <- gsub("antibody", "", colnames(designlimmaP))

contrastP <- limma::makeContrasts(T5vsT0 = (PD1_5-Con_5)-(PD1_0-Con_0),
                                  T20vsT0 = (PD1_20-Con_20)-(PD1_0-Con_0),
                                  T24vsT0 = (PD1_24-Con_24)-(PD1_0-Con_0),
                                  T20vsT5 = (PD1_20-Con_20)-(PD1_5-Con_5),
                                  T24vsT5 = (PD1_24-Con_24)-(PD1_5-Con_5),
                                  T24vsT20 = (PD1_24-Con_24)-(PD1_20-Con_20),
                                  PD1_0vsCon0 = PD1_0 - Con_0,
                                  PD1_0vsCon5 = PD1_0 - Con_5,
                                  PD1_0vsCon20 = PD1_0 - Con_20,
                                  PD1_0vsCon24 = PD1_0 - Con_24,
                                  PD1_5vsCon0 = PD1_5 - Con_0,
                                  PD1_5vsCon5 = PD1_5 - Con_5,
                                  PD1_5vsCon20 = PD1_5 - Con_20,
                                  PD1_5vsCon24 = PD1_5 - Con_24,
                                  PD1_20vsCon0 = PD1_20 - Con_0,
                                  PD1_20vsCon5 = PD1_20 - Con_5,
                                  PD1_20vsCon20 = PD1_20 - Con_20,
                                  PD1_20vsCon24 = PD1_24 - Con_24,
                                  PD1_24vsCon0 = PD1_24 - Con_0,
                                  PD1_24vsCon5 = PD1_24 - Con_5,
                                  PD1_24vsCon20 = PD1_24 - Con_20,
                                  PD1_24vsCon24 = PD1_24 - Con_24,
                                  levels=colnames(designlimmaP))

## fitting model ####
fit_interaction <- lmFit(imputed_nostim, designlimmaP) # to check
fit_interaction <- contrasts.fit(fit_interaction, contrastP)
fit_interaction <- eBayes(fit_interaction) 

coef_ctrl_ant <- list() ### varian of code: https://www.r-bloggers.com/concatenating-a-list-of-data-frames/
for (s in 1:length(colnames(contrastP))) {
  coef_ctrl_ant[[s]] <- topTable(fit_interaction, coef=s, n=Inf)
  setDT(coef_ctrl_ant[[s]], keep.rownames = TRUE)[]
  coef_ctrl_ant[[s]]$type <- colnames(contrastP)[s]
}

diff_exp_protein <- do.call(rbind, coef_ctrl_ant)

## filter based on adj pval (<= 0.05) ####
diff_exp_protein <- diff_exp_protein[which(diff_exp_protein$adj.P.Val <= .05),]
## filter based on LogFC (>= 2) ####
diff_exp_protein <- diff_exp_protein[which(diff_exp_protein$logFC >= 2),]

## adding gene names to DEP ####
gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% diff_exp_protein$rn),]
diff_exp_protein <- diff_exp_protein %>%
  left_join(gene_name_matrix[,c(1,4)], by = c("rn" = "Protein.IDs"))

## DEP by comparison type #### 
time_genes <- diff_exp_protein[which(diff_exp_protein$type %in% c("T5vsT0", "T20vsT0", "T24vsT0", 
                                                                  "T20vsT5",  "T24vsT5", "T24vsT20")),]
                                     
interaction_genes <- diff_exp_protein[which(diff_exp_protein$type %in% c("PD1_0vsCon0", "PD1_0vsCon5",
                                                                         "PD1_0vsCon20","PD1_0vsCon24", 
                                                                         "PD1_5vsCon0", "PD1_5vsCon5",
                                                                         "PD1_5vsCon20", "PD1_5vsCon24", 
                                                                         "PD1_20vsCon0", "PD1_20vsCon5", 
                                                                         "PD1_20vsCon20", "PD1_20vsCon24", 
                                                                         "PD1_24vsCon0", "PD1_24vsCon5", 
                                                                         "PD1_24vsCon20", "PD1_24vsCon24")),]


gene_list <- list(time = time_genes$Gene.names,
                  interaction = interaction_genes$Gene.names)

gene_list <- Venn(gene_list)
time_genes_only <- discern(gene_list, "time", "interaction")
interaction_genes_only <- discern(gene_list, "interaction", "time")
overlap_int_time <- overlap(gene_list)

time_genes_only.df <- time_genes[which(time_genes$Gene.names %in% time_genes_only),]
time_genes_only.df <- time_genes_only.df[!duplicated(time_genes_only.df[,"Gene.names"]),]

interaction_genes_only.df <- interaction_genes[which(interaction_genes$Gene.names %in% interaction_genes_only),]
interaction_genes_only.df <- interaction_genes_only.df[!duplicated(interaction_genes_only.df[,"Gene.names"]),]

ovelapped_genes.df <- diff_exp_protein[which(diff_exp_protein$Gene.names %in% overlap_int_time),]
ovelapped_genes.df <- ovelapped_genes.df[!duplicated(ovelapped_genes.df[,"Gene.names"]),]


ggvenn(gene_list)   

## clustering time genes ####
imputed_filtered_time <- imputed_nostim[which(rownames(imputed_nostim) %in% time_genes_only.df$rn),]

gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(imputed_filtered_time)),]
imputed_filtered_time  <- cbind(imputed_filtered_time , gene_name_matrix[, "Gene.names", drop=FALSE])
imputed_filtered_time  <- imputed_filtered_time [!duplicated(imputed_filtered_time[,"Gene.names"]),]
rownames(imputed_filtered_time) <- NULL
imputed_filtered_time <- column_to_rownames(imputed_filtered_time , "Gene.names")

annotation_nostim$condition <- substr(annotation_nostim$antibody, 1, 3)

pdf("/home/degan/ip_proteomics/figures/interaction_term/clustering_time.pdf", height = , width = 6)
pheatmap(imputed_filtered_time, cluster_rows = T, show_colnames = F,
         fontsize_row = 5, scale = "row", annotation = annotation_nostim[,c(1,2,3),])
dev.off()

## clustering antibody genes ####
imputed_filtered_time <- imputed_nostim[which(rownames(imputed_nostim) %in% antibody_genes$rn),]

gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(imputed_filtered_time)),]
imputed_filtered_time  <- cbind(imputed_filtered_time , gene_name_matrix[, "Gene.names", drop=FALSE])
rownames(imputed_filtered_time) <- NULL
imputed_filtered_time <- column_to_rownames(imputed_filtered_time , "Gene.names")

pheatmap(imputed_filtered_time, cluster_rows = T, show_colnames = F,
         fontsize_row = 6, scale = "row", annotation = annotation_nostim[,2, drop=F])

## clustering interaction genes ####
imputed_filtered_time <- imputed_nostim[which(rownames(imputed_nostim) %in% interaction_genes$rn),]

gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(imputed_filtered_time)),]
imputed_filtered_time  <- cbind(imputed_filtered_time , gene_name_matrix[, "Gene.names", drop=FALSE])
imputed_filtered_time  <- imputed_filtered_time [!duplicated(imputed_filtered_time[,"Gene.names"]),]
rownames(imputed_filtered_time) <- NULL
imputed_filtered_time <- column_to_rownames(imputed_filtered_time , "Gene.names")

pheatmap(imputed_filtered_time, cluster_rows = T, show_colnames = F,
         fontsize_row = 6, scale = "row", annotation = annotation_nostim)


## fgsea ####

ensembl.pathway <- sbgn.gsets(id.type = "SYMBOL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = T, #T
                              #database = "MetaCyc", 
                              truncate.name.length = 100)

genes <- ovelapped_genes.df %>% pull(logFC, Gene.names)

fgseaRes <- fgsea(pathways = ensembl.pathway, stats = genes)
fgseaRes$pathway <- gsub(".*::","",fgseaRes$pathway)
fgseaRes <- fgseaRes %>% arrange(desc(abs(NES)))



