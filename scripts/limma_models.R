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

## read in exp design - created by Vadim ####
experimental_design <- read.csv(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/experimental_design1.csv", sep = "\t", 
                                colClasses=c("character","character"), header = F)
colnames(experimental_design) <- c("label","condition")

##remove stimulation from exp design
experimental_design$stimulation <- ifelse(grepl("Stim",experimental_design$label), "Stimulated", "NotStim")
experimental_design <- experimental_design[which(experimental_design$stimulation == "NotStim"),]


## creating design matrix ####
groups <- as.factor(experimental_design$condition)
designlimmaP <- as.data.frame(model.matrix(~ 0 + groups))
names(designlimmaP) <- gsub("groups", "", names(designlimmaP))

## contrast each control and PD1 for each antibody and timepoint  ####  
contrastP <- limma::makeContrasts(paste(levels(groups)[1],"-",levels(groups)[13],sep = ""),
                    paste(levels(groups)[2],"-",levels(groups)[14],sep = ""),
                    paste(levels(groups)[3],"-",levels(groups)[15],sep = ""),
                    paste(levels(groups)[4],"-",levels(groups)[16],sep = ""),
                    
                    
                    paste(levels(groups)[5],"-",levels(groups)[17],sep = ""),
                    paste(levels(groups)[6],"-",levels(groups)[18],sep = ""),
                    paste(levels(groups)[7],"-",levels(groups)[19],sep = ""),
                    paste(levels(groups)[8],"-",levels(groups)[20],sep = ""),
                    paste(levels(groups)[9],"-",levels(groups)[21],sep = ""),
                    
                    paste(levels(groups)[10],"-",levels(groups)[22],sep = ""),
                    paste(levels(groups)[11],"-",levels(groups)[23],sep = ""),
                    paste(levels(groups)[12],"-",levels(groups)[24],sep = ""),
                    levels=designlimmaP)

fit_ctrl_ant <- lmFit(imputed_nostim, designlimmaP) # to check
fit_ctrl_ant <- contrasts.fit(fit_ctrl_ant, contrastP)
fit_ctrl_ant<- eBayes(fit_ctrl_ant) 

coef_ctrl_ant <- list() ### varian of code: https://www.r-bloggers.com/concatenating-a-list-of-data-frames/
for (s in 1:length(colnames(contrastP))) {
  coef_ctrl_ant[[s]] <- topTable(fit_ctrl_ant, coef=s, n=Inf)
  setDT(coef_ctrl_ant[[s]], keep.rownames = TRUE)[]
  coef_ctrl_ant[[s]]$type <- colnames(contrastP)[s]
}

diff_exp_protein <- do.call(rbind, coef_ctrl_ant)

## filter based on adj pval (<= 0.05) ####
diff_exp_protein <- diff_exp_protein[which(diff_exp_protein$adj.P.Val <= .05),]
## filter based on LogFC (>= 2) ####
diff_exp_protein <- diff_exp_protein[which(diff_exp_protein$logFC >= 2),]

## adding gene names to d_e_p ####
gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% diff_exp_protein$rn),]
diff_exp_protein <- diff_exp_protein %>%
  left_join(gene_name_matrix[,c(1:4)], by = c("rn" = "Protein.IDs"))

## number of diff exp proteins and plotting as barplot ####
number_proteins <- data.frame(table(diff_exp_protein$type))
number_proteins$antibody <- sapply(strsplit(as.character(number_proteins$Var1), "_"), "[", 1)
number_proteins$antibody <- str_sub(number_proteins$antibody, 4, length(number_proteins$antibody))

pdf("/home/degan/ip_proteomics/figures/boxplot_freq_diff_exp.pdf", height = 5, width = 4)
ggplot(data=number_proteins, aes(x=Var1, y=Freq, fill = antibody)) + 
  geom_bar(stat="identity") + scale_fill_brewer(type = "qual", palette = "Set2") + theme_bw() + xlab("Condition comparison") + ylab("Differentially Exp Proteins (Freq)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
dev.off()

## the number conditions a protein is differentially exp in ####
## remove proteins diff exp in only 1 condition and carry out gsea ####
protein_numbers <- data.frame(table(diff_exp_protein$Gene.names))
protein_numbers <- protein_numbers[which(protein_numbers$Freq > 1),]
diff_exp_protein_filtered <- diff_exp_protein[which(diff_exp_protein$Gene.names %in% protein_numbers$Var1),]

ensembl.pathway <- sbgn.gsets(id.type = "SYMBOL",
                              species = "hsa",
                              mol.type = "gene",
                              output.pathway.name = T, #T
                              #database = "MetaCyc", 
                              truncate.name.length = 100)

genes <- diff_exp_protein_filtered %>% pull(logFC, Gene.names)

fgseaRes <- fgsea(pathways = ensembl.pathway, stats = genes)
fgseaRes$pathway <- gsub(".*::","",fgseaRes$pathway)
fgseaRes <- fgseaRes %>% arrange(desc(abs(NES)))

## creating upset plot for each condition ####
col_bars <- brewer.pal(n=3, name = "Set1")
upset_df <- diff_exp_protein
upset_df$bin <- ifelse(diff_exp_protein$logFC >=2 & diff_exp_protein$adj.P.Val <0.05, 1, 0)
upset_df <- acast(upset_df, type~rn, value.var="bin")
upset_df <- t(upset_df)
upset_df <- rownames_to_column(data.frame(upset_df))
names(upset_df)[1] <- "Name"
pdf("/home/degan/ip_proteomics/figures/upsetplot.pdf", height = 5, width = 7)
upset(upset_df, sets = colnames(upset_df)[2:13], order.by = "freq", 
      nsets = 12, keep.order = T, point.size = 1.8, text.scale = .8)
dev.off() 

## adding time as a model variable ####
condition <- substr(sapply(strsplit(colnames(imputed_nostim), "_"), "[", 1), 1, 3)
antibody <- substr(sapply(strsplit(colnames(imputed_nostim), "_"), "[", 1), 4,6)
time <-  sapply(strsplit(colnames(imputed_nostim), "_"), "[", 2)

time <- gsub("min", "", time)
time <- gsub("hr", "", time)
time <- ifelse(time == 24, as.numeric(time)*60, as.numeric(time)*1)

condition_time <- paste(condition,time, sep="_")

designlimmaP <- as.data.frame(model.matrix(~ 0 + condition_time))
colnames(designlimmaP) <- gsub("condition_time", "", colnames(designlimmaP))

cor <- duplicateCorrelation(imputed_nostim, designlimmaP, block=antibody)
cor$consensus.correlation

fit <- lmFit(object=imputed_nostim, design=designlimmaP, 
             block=antibody, correlation=cor$consensus.correlation)
fit$coefficients[1:5,]

## contrast each control and PD1 for each antibody and timepoint  ####  
contrastP <- limma::makeContrasts(PD1_0vsCON_0 = PD1_0 - Con_0,
                                  PD1_5vsCON_5 = PD1_5 - Con_5,
                                  PD1_20vsCON_20 = PD1_20 - Con_20,
                                  PD1_24vsCON_24 = PD1_24 - Con_24,
                                  levels=colnames(designlimmaP))

fit_ant_blocked <- contrasts.fit(fit, contrastP)
fit_ant_blocked <- eBayes(fit_ant_blocked) 
tt.PD1vsCtr <- topTable(fit_ant_blocked, number = nrow(fit_ant_blocked))

coef_ctrl_ant <- list() ### varian of code: https://www.r-bloggers.com/concatenating-a-list-of-data-frames/
for (s in 1:length(colnames(contrastP))) {
  coef_ctrl_ant[[s]] <- topTable(fit_ant_blocked, coef=s, n=Inf)
  setDT(coef_ctrl_ant[[s]], keep.rownames = TRUE)[]
  coef_ctrl_ant[[s]]$type <- colnames(contrastP)[s]
}

diff_exp_protein <- do.call(rbind, coef_ctrl_ant)
## filter based on adj pval (<= 0.05) ####
diff_exp_protein <- diff_exp_protein[which(diff_exp_protein$adj.P.Val <= .05),]
## filter based on LogFC (>= 2) ####
diff_exp_protein <- diff_exp_protein[which(diff_exp_protein$logFC >= 1.5),]

## adding gene names to d_e_p ####
gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% diff_exp_protein$rn),]
diff_exp_protein <- diff_exp_protein %>%
  left_join(gene_name_matrix[,c("Protein.IDs","Gene.names")], by = c("rn" = "Protein.IDs"))

## remove PDCD1; we know its pulled down, and throws off significane level in plots ####
diff_exp_protein <- diff_exp_protein[which(!diff_exp_protein$Gene.names == "PDCD1"),]

## plotting differentially expressed proteins
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

## plotting DEP in grid format (volcano plot) #### 
library("gridExtra")
pdf("/home/degan/ip_proteomics/figures/volcano_plot_PD1vsctr.pdf", height = 4, width = 10)
grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
             plot_list[[4]], nrow =1)
dev.off()

## saving each volcano plot individually ####
pdf("/home/degan/ip_proteomics/figures/volcano_plot_PD1vsctr0.pdf", height = 4, width = 6)
plot_list[[1]]
dev.off()

pdf("/home/degan/ip_proteomics/figures/volcano_plot_PD1vsctr2.pdf", height = 4, width = 6)
plot_list[[2]]
dev.off()

pdf("/home/degan/ip_proteomics/figures/volcano_plot_PD1vsctr3.pdf", height = 4, width = 6)
plot_list[[3]]
dev.off()

pdf("/home/degan/ip_proteomics/figures/volcano_plot_PD1vsctr4.pdf", height = 4, width = 6)
plot_list[[4]]
dev.off()

## plotting upset plot for PD1vsCtr ####
upset_df <- diff_exp_protein
upset_df$bin <- ifelse(diff_exp_protein$logFC >=1.5 & diff_exp_protein$adj.P.Val <0.05, 1, 0)
upset_df <- acast(upset_df, type~rn, value.var="bin")
upset_df <- t(upset_df)
upset_df <- rownames_to_column(data.frame(upset_df))
names(upset_df)[1] <- "Name"
pdf("/home/degan/ip_proteomics/figures/upsetplot_pd1_vs_ctr.pdf", height = 5, width = 7)
upset(upset_df, sets = colnames(upset_df)[2:5], order.by = "freq",
      keep.order = T, point.size = 1.8, text.scale = .8)
dev.off()



## adding gene names to d_e_p ####
gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(tt.PD1vsCtr)),]
tt.PD1vsCtr <- tt.PD1vsCtr %>% merge(gene_name_matrix[,c(1:4)], by.x = 0, by.y = "Protein.IDs")

tt.PD1vsCtr_filtered <- tt.PD1vsCtr[which(tt.PD1vsCtr$adj.P.Val < 0.005),]
tt.PD1vsCtr_filtered  <- tt.PD1vsCtr_filtered %>% arrange(desc(abs(F)))
tt.PD1vsCtr_filtered <- tt.PD1vsCtr_filtered[1:100, ]

## Building gene modules using filtered genes ####
set.seed(16351893)
imputed_filtered_protein <- imputed_nostim[which(rownames(imputed_nostim) %in% tt.PD1vsCtr_filtered$Row.names),]

heat_annotation <- ComplexHeatmap::HeatmapAnnotation(experiment = substr(annotation_nostim$antibody, 1, 3),
                                                     condition = annotation_nostim$antibody,
                                                     time = annotation_nostim$time,
                                                     show_annotation_name = FALSE)

cor_matrix <- cor(imputed_filtered_protein)
pdf("/home/degan/ip_proteomics/figures/heatmap_antibodytype_blocked.pdf", width = 5, height = 5)
ComplexHeatmap::Heatmap(cor_matrix,name = "Corr", top_annotation = heat_annotation, show_row_names = FALSE, 
                        show_column_names = FALSE)
saveRDS(imputed_filtered_protein, "/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")
dev.off()




