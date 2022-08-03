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
library(openxlsx)
library(tidyverse)

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

## saving differentially bound proteins as xls file ####
write.xlsx(diff_exp_protein, '/home/degan/ip_proteomics/inputs/differentially_bound_proteins_long.xlsx')
saveRDS(diff_exp_protein, '/home/degan/ip_proteomics/inputs/differentially_bound_proteins.Rds')
## creating upset plot for each condition ####
col_bars <- brewer.pal(n=3, name = "Set1")
upset_df <- diff_exp_protein
upset_df$bin <- ifelse(diff_exp_protein$logFC > 1.5 & diff_exp_protein$adj.P.Val <0.05, 1, 0)
upset_df <- acast(upset_df, type~rn, value.var="bin")
upset_df <- t(upset_df)
upset_df <- rownames_to_column(data.frame(upset_df))
names(upset_df)[1] <- "Name"
pdf("/home/degan/ip_proteomics/figures/antibody_fixed/upsetplot.pdf", height = 5, width = 7)
upset(upset_df, sets = colnames(upset_df)[2:9], order.by = "freq", 
      nsets = 12, keep.order = T, point.size = 1.8, text.scale = .8)
dev.off() 


## proteins diff expressed all time points  ####
common_genes <- diff_exp_protein[diff_exp_protein$type %in% c("PD1_0vsCon_0", "PD1_20vsCon_20", "PD1_24vsCon_24",   "PD1_5vsCon_5"),]
common_genes$bin <- ifelse(abs(common_genes$logFC) > 1 & common_genes$adj.P.Val <0.05, 1, 0)
common_genes <- acast(common_genes, type~rn, value.var="bin")
common_genes <- t(common_genes)
common_genes <- rownames_to_column(data.frame(common_genes))
common_genes$sum <- rowSums(common_genes[,c(2:length(colnames(common_genes)))])
common_genes<- common_genes %>%
  left_join(gene_name_matrix[,c("Protein.IDs", "Gene.names")], by = c("rowname" = "Protein.IDs"))
common_genes <- common_genes[common_genes$sum == 4,]

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
                       labSize = 2,
                       pointSize = 1,
                       max.overlaps = 30,
                       #transcriptLabSize = 4,
                       #transcriptLabhjust = 0.5,
                       #legend=c("","","",""),
                       #legendPosition = "right",
                       legendPosition = 'none'
  )
  t = max(temp$logFC)
  plot_list[[i]] = p + ggplot2::coord_cartesian(xlim=c(-t, t)) 
  
}

#plotting time-course 
## saving each volcano plot individually ####
pdf("/home/degan/ip_proteomics/figures/antibody_fixed/volcano_timecourse_pd1/volcano_plot_PD1vsCtr0.pdf", height = 5, width = 6)
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

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/volcano_timecourse_pd1/volcano_plot_PD1vsCtr.pdf", height = 5, width = 6)
plot_list[["PD1vsCtr"]]
dev.off()

## Filter proteins according to logFC and p-value ####
diff_exp_filter <- diff_exp_protein[which(abs(diff_exp_protein$logFC) > 1 & diff_exp_protein$adj.P.Val <0.05),]
saveRDS(diff_exp_filter, "/home/degan/ip_proteomics/inputs/diff_bound_filter.Rds")

## save as wide format ####
diff_exp_wide <- diff_exp_filter[,c("logFC", "Gene.names", "type")]
diff_exp_wide <- reshape(diff_exp_wide, idvar = "Gene.names", timevar = "type", direction = "wide")
diff_exp_wide[is.na(diff_exp_wide)] <- 0
rownames(diff_exp_wide) <- NULL
diff_exp_wide <- column_to_rownames(diff_exp_wide, "Gene.names")

write.xlsx(diff_exp_wide, '/home/degan/ip_proteomics/inputs/differentially_bound_proteins_filter_wide.xlsx')

#gene_list <- list(PD1vsCtr = diff_exp_filter$Gene.names[diff_exp_filter$type=="PD1vsCtr"],
#                  "PD1-time" = diff_exp_filter$Gene.names[diff_exp_filter$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")],
#                  "PD1vsCtr-time" = diff_exp_filter$Gene.names[diff_exp_filter$type %in% c("PD1_0vsCon_0", "PD1_20vsCon_20", 
#                                                                       "PD1_24vsCon_24",   "PD1_5vsCon_5")])

#gene_list <- data.frame(comparison = diff_exp_filter$type, gene=diff_exp_filter$Gene.names)
#gene_list <- gene_list  %>% group_by(comparison) %>% summarise(across(everything(), str_c, collapse=" ")) 
#gene_list <- gene_list[c(1:4),]
#gene_list$comparison <- factor(gene_list$comparison,levels=c("PD1_0vsCon_0", "PD1_5vsCon_5","PD1_20vsCon_20", "PD1_24vsCon_24"))
#gene_list$gene <- strsplit(gene_list$gene, split = " ")
#gene_list_1 <- gene_list$gene
#names(gene_list_1) <- gene_list$comparison

gene_list_v <- Venn(gene_list)
pdf("/home/degan/ip_proteomics/figures/antibody_fixed/venn_diagram.pdf")
ggvenn(gene_list_v)
dev.off()

## fgsea for each comparison ###
msigdbr_df <- msigdbr(species = "human", category = "C2")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

ensembl.pathway <- sbgn.gsets(id.type = "SYMBOL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = T, #T
                              database = "MetaCyc", 
                              truncate.name.length = 100)

top_proteins_pd1time <- diff_exp_filter[(diff_exp_filter$Gene.names %in% gene_list[["PD1-time"]] &
                                         diff_exp_filter$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")),] %>% 
                                         arrange(desc(t))

top_proteins_pd1_vsctr <- diff_exp_filter[(diff_exp_filter$Gene.names %in% gene_list[["PD1vsCtr-time"]] &
                                           diff_exp_filter$type %in% c("PD1_24vsCon_24")),] %>% 
  arrange(desc(t))
saveRDS(top_proteins_pd1time, "/home/degan/ip_proteomics/inputs/ip_proteins_pd1time.Rds")

genes <- top_proteins_pd1_vsctr %>% pull(logFC, Gene.names)
names(genes) <- sub(";.*", "", names(genes))

fgseaRes <- fgsea(pathways = ensembl.pathway, stats = genes, nPermSimple = 10000)


## number of diff exp proteins and plotting as barplot ####
number_proteins <- data.frame(table(diff_exp_filter$type))
number_proteins <- number_proteins[c(1:4),]
number_proteins$Var1 <- factor(number_proteins$Var1,levels=c("PD1_0vsCon_0", "PD1_5vsCon_5","PD1_20vsCon_20", "PD1_24vsCon_24"))

pdf("/home/degan/ip_proteomics/figures/antibody_fixed/boxplot_freq_diff_exp_time.pdf", height = 3.5, width = 4)
ggplot(data=number_proteins, aes(x=Var1, y=Freq, fill = Var1)) + 
  geom_bar(stat="identity") + scale_fill_npg() + theme_bw() + xlab("Comparison") + ylab("Diff Bound Proteins (Log2FC>1, adj.p.val<0.05)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10)) + labs(fill="") 
dev.off()


## creating upset plot for each time point ####
upset_df <- diff_exp_protein
upset_df <- upset_df[diff_exp_protein$type %in% c("PD1_0vsCon_0", "PD1_5vsCon_5","PD1_20vsCon_20", "PD1_24vsCon_24"),]
upset_df$type <- factor(upset_df$type,levels=c("PD1_0vsCon_0", "PD1_5vsCon_5","PD1_20vsCon_20", "PD1_24vsCon_24"))
upset_df$bin <- ifelse(abs(upset_df$logFC) > 1 & upset_df$adj.P.Val <0.05, 1, 0)
upset_df <- acast(upset_df, type~rn, value.var="bin")
upset_df <- t(upset_df)
upset_df <- rownames_to_column(data.frame(upset_df))
names(upset_df)[1] <- "Name"
pdf("/home/degan/ip_proteomics/figures/antibody_fixed/upsetplot_timepoints.pdf", height = 3, width = 4)
upset(upset_df, sets = colnames(upset_df)[2:5], order.by = "freq",
      nsets = 12, keep.order = T, point.size = 1.8, text.scale = .8, matrix.color = pal_npg("nrc")(10)[2])
dev.off() 
























