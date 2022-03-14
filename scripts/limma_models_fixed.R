################################################################################
# Donagh Egan
# POIAZ 
# Date: 2nd March 2022
# 

# Description: This script runs a linear model comparing conditions for a PD1 
# immunoprecipitation study in a Jurkat cell line.         
################################################################################

## Library #### 
################################################################################

library(limma)
library(data.table)
library(dplyr)
library(stringr)
library(fgsea)
library(msigdbr)
library(SBGNview)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(reshape2)
library(RColorBrewer)
library(stringi)
library(UpSetR)
library(RVenn)

## Design 1: ####
## Comparing PD1 pull down exps with cd3/cd28 stimulation experiments across time points ####

## Read in imputed data matrix and anntotations ####
pgroups <- read.table(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/proteinGroups_PPI_Martina_Jan.txt", 
                      header = T, sep = "\t",quote='',stringsAsFactors = FALSE,comment.char="")
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

## model variables ####
condition <- annotation_df$condition
antibody <- annotation_df$antibody
time <- annotation_df$time
stim <- annotation_df$stimulation

## formatting variables ####
time <- stri_replace_all_regex(time, pattern = c('min', 'hr'), replacement = "", vectorize = F)
condition_time <- paste(condition,time, sep="_") ## combining time and condition 
con_time_stim <- paste(condition_time, stim, sep = "_")

## design matrix ####
designlimmaP <- as.data.frame(model.matrix(~ 0 + con_time_stim + antibody))
names(designlimmaP) <- gsub("con_time_stim", "", names(designlimmaP))

## contrasts ####
# proteins that are differentially bound as a function of time  
contrastP <- limma::makeContrasts(PD1vsCtr = (PD1_0_unstim + PD1_5_stim + 
                                              PD1_20_stim + PD1_24_stim)/4 - 
                                             (Con_0_unstim + Con_5_stim +
                                              Con_20_stim + Con_24_stim)/4,
                                  PD1_0vsCon_0 = PD1_0_unstim - Con_0_unstim,
                                  PD1_5vsCon_5 = PD1_5_stim-Con_5_stim,
                                  PD1_20vsCon_20 = PD1_20_stim-Con_20_stim,
                                  PD1_24vsCon_24 = PD1_24_stim-Con_24_stim,
                                  ST5vNsT0 = (PD1_5_stim-Con_5_stim)-(PD1_0_unstim-Con_0_unstim),
                                  ST20vST5 = (PD1_20_stim-Con_20_stim) - (PD1_5_stim-Con_5_stim),
                                  ST24vST20 = ((PD1_24_stim-PD1_24_unstim)-(Con_24_stim-Con_0_unstim)) - (PD1_20_stim-Con_20_stim),
                                  levels=colnames(designlimmaP))

## fitting model to data and making contrasts ####
fit_ctrl_ant <- lmFit(imputed_data, designlimmaP) # to check
fit_ctrl_ant <- contrasts.fit(fit_ctrl_ant, contrastP)
fit_ctrl_ant<- eBayes(fit_ctrl_ant) 

## combining differentially expressed proteins ####
coef_ctrl_ant <- list() ### variant of code: https://www.r-bloggers.com/concatenating-a-list-of-data-frames/
for (s in 1:length(colnames(contrastP))) {
  coef_ctrl_ant[[s]] <- topTable(fit_ctrl_ant, coef=s, n=Inf)
  setDT(coef_ctrl_ant[[s]], keep.rownames = TRUE)[]
  coef_ctrl_ant[[s]]$type <- colnames(contrastP)[s]
}

diff_exp_protein <- do.call(rbind, coef_ctrl_ant)

## adding gene names to d_e_p ####
gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% diff_exp_protein$rn),]
diff_exp_protein <- diff_exp_protein %>%
  left_join(gene_name_matrix[,c("Protein.IDs", "Gene.names")], by = c("rn" = "Protein.IDs"))

## creating upset plot for each condition ####
col_bars <- brewer.pal(n=3, name = "Set1")
upset_df <- diff_exp_protein
upset_df$bin <- ifelse(diff_exp_protein$logFC >=1.5 & diff_exp_protein$adj.P.Val <0.05, 1, 0)
upset_df <- acast(upset_df, type~rn, value.var="bin")
upset_df <- t(upset_df)
upset_df <- rownames_to_column(data.frame(upset_df))
names(upset_df)[1] <- "Name"
pdf("/home/degan/ip_proteomics/figures/antibody_fixed/upsetplot.pdf", height = 5, width = 7)
upset(upset_df, sets = colnames(upset_df)[2:9], order.by = "freq", 
      nsets = 12, keep.order = T, point.size = 1.8, text.scale = .8)
dev.off() 


## proteins diff expressed all time points  ####
common_genes <- diff_exp_protein
common_genes$bin <- ifelse(common_genes$logFC >=1.5 & common_genes$adj.P.Val <0.05, 1, 0)
common_genes <- acast(common_genes, type~rn, value.var="bin")
common_genes <- t(common_genes)
common_genes <- rownames_to_column(data.frame(common_genes))
common_genes$sum <- rowSums(common_genes[,c(2:length(colnames(common_genes)))])

## plotting differentially expressed proteins ####
library(EnhancedVolcano)
g1 <- levels(as.factor(diff_exp_protein$type))
plot_list = list()
for (i in g1) {
  temp <- subset(diff_exp_protein, type == paste0(i))
  temp <- as.data.frame(temp[,c(2,6,8,9)])
  p <- EnhancedVolcano(temp,
                       lab = temp$Gene.names,
                       x = "logFC",
                       y = "adj.P.Val",
                       pCutoff = 0.05, 
                       FCcutoff = log2(1.5),
                       col=c("grey", "grey", "grey", "red3"),
                       colAlpha = 1, subtitle = " ",
                       titleLabSize = 12,
                       axisLabSize = 12,
                       title = paste(i),
                       captionLabSize = 0,
                       labSize = 2.5,
                       pointSize = 1,
                       max.overlaps = 20,
                       #transcriptLabSize = 4,
                       #transcriptLabhjust = 0.5,
                       #legend=c("","","",""),
                       #legendPosition = "right",
                       legendPosition = 'none'
  )
  t <- max(temp$logFC) + 1
  
  plot_list[[i]] = p + ggplot2::coord_cartesian(xlim=c(-t, t)) 
  
}

#plotting time-course 
## saving each volcano plot individually ####
pdf("/home/degan/ip_proteomics/figures/antibody_fixed/volcano_timecourse_pd1/volcano_plot_PD1vsCtr.pdf", height = 5, width = 6)
plot_list[[1]]
dev.off()

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/volcano_timecourse_pd1/volcano_plot_ST20vST5.pdf", height = 5, width = 6)
plot_list[[2]]
dev.off()

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/volcano_timecourse_pd1/volcano_plot_ST24vST20.pdf", height = 5, width = 6)
plot_list[[3]]
dev.off()

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/volcano_timecourse_pd1/volcano_plot_ST5vNsT0.pdf", height = 5, width = 6)
plot_list[[7]]
dev.off()

## clustering based on top proteins ####
## top 100 unique d_e_p according to t value  ####
set.seed(16351893)
top_proteins <- diff_exp_protein %>% arrange(desc(t)) %>% group_by(type) %>% 
              dplyr::slice(1:10)

annotation_subset <- annotation_df[which(annotation_df$time == "24hr" & annotation_df$stimulation == "unstim" |
                                           annotation_df$condition == "Con"),]
annotation_subset <- annotation_df[which(!rownames(annotation_df) %in% rownames(annotation_subset)),]

imputed_top <- imputed_data[which(rownames(imputed_data) %in% top_proteins$rn),]


imputed_top <- imputed_top[,which(colnames(imputed_top) %in% rownames(annotation_subset))]

## resolving gene names ####
imputed_top <- merge(top_proteins[,c("rn", "Gene.names")], imputed_top, by.x = "rn", by.y=0)
imputed_top <- imputed_top[!duplicated(imputed_top[,"Gene.names"]),]
rownames(imputed_top) <- NULL 
imputed_top <- imputed_top %>% column_to_rownames("Gene.names")
imputed_top$rn <- NULL

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/timecourse_PD1_clustering.pdf")
pheatmap::pheatmap(imputed_top, scale = "row", show_colnames = FALSE, show_rownames = T,
                   annotation = annotation_df[,c(1,3,5), drop=F], treeheight_row = 15,
                   treeheight_col = 25, fontsize_row = 6)
dev.off()


## venn diagram ####
diff_exp_filter <- diff_exp_protein[which(diff_exp_protein$logFC >=1.5 & diff_exp_protein$adj.P.Val <0.05),]
saveRDS(diff_exp_filter, "/home/degan/ip_proteomics/inputs/diff_bound_filter.Rds")

gene_list <- list(PD1vsCtr = diff_exp_filter$Gene.names[diff_exp_filter$type=="PD1vsCtr"],
                  "PD1-time" = diff_exp_filter$Gene.names[diff_exp_filter$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")],
                  "PD1vsCtr-time" = diff_exp_filter$Gene.names[diff_exp_filter$type %in% c("PD1_0vsCon_0", "PD1_20vsCon_20", 
                                                                       "PD1_24vsCon_24",   "PD1_5vsCon_5")])
gene_list_v <- Venn(gene_list)
pdf("/home/degan/ip_proteomics/figures/antibody_fixed/venn_diagram.pdf")
ggvenn(gene_list_v)
dev.off()

## clustering according to pd1 time-course genes ####
set.seed(16351893)
top_proteins_pd1time <- diff_exp_filter[(diff_exp_filter$Gene.names %in% gene_list[["PD1-time"]] &
                                         diff_exp_filter$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")),] %>% arrange(desc(t)) %>% 
                                         group_by(type) %>% dplyr::slice(1:10)

annotation_subset <- annotation_df[which(!annotation_df$condition == "Con"),]
annotation_subset <- annotation_subset[!(annotation_subset$stimulation == "unstim" &
                                         annotation_subset$time == "24hr"),]

imputed_top_pd1time <- imputed_data[which(rownames(imputed_data) %in% top_proteins_pd1time$rn),]


## resolving gene names ####
imputed_top_pd1time <- merge(top_proteins_pd1time[,c("rn", "Gene.names")], imputed_top_pd1time, by.x = "rn", by.y=0)
imputed_top_pd1time <- imputed_top_pd1time[!duplicated(imputed_top_pd1time[,"Gene.names"]),]
rownames(imputed_top_pd1time) <- NULL 
imputed_top_pd1time <- imputed_top_pd1time %>% column_to_rownames("Gene.names")
imputed_top_pd1time$rn <- NULL

imputed_top_pd1time <- imputed_top_pd1time[,which(colnames(imputed_top_pd1time) %in% rownames(annotation_subset))]

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/timecourse_PD1_clustering.pdf")
pheatmap::pheatmap(imputed_top_pd1time, scale = "row", show_colnames = FALSE, show_rownames = T,
                   annotation = annotation_subset[,c(1,3,5), drop=F], treeheight_row = 15,
                   treeheight_col = 25, fontsize_row = 6)
dev.off()

## fgsea time-course ###
msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

ensembl.pathway <- sbgn.gsets(id.type = "SYMBOL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = T, #T
                              #database = "MetaCyc", 
                              truncate.name.length = 100)

top_proteins_pd1time <- diff_exp_filter[(diff_exp_filter$Gene.names %in% gene_list[["PD1-time"]] &
                                         diff_exp_filter$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")),] %>% 
                                         arrange(desc(t))
saveRDS(top_proteins_pd1time, "/home/degan/ip_proteomics/inputs/ip_proteins_pd1time.Rds")

genes <- top_proteins_pd1time %>% pull(logFC, Gene.names)
names(genes) <- sub(";.*", "", names(genes))

fgseaRes <- fgsea(pathways = pathwaysH, stats = genes)






















