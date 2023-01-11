library(dplyr)
library(tidyverse)
library(ggsci)
library(fgsea)
library(msigdbr)
library(ggpubr)

cluster_info <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
pd1_proteins <- cluster_info[cluster_info$cluster_id == "PD1/NF-kappaB",]
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
gene_names <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds"))

pd1_proteins <- imputed_data[rownames(imputed_data) %in% pd1_proteins$Row.names,]
pd1_proteins <- t(pd1_proteins)
pd1_proteins <- merge(pd1_proteins, annotation_df, by =0)
pd1_proteins <- pd1_proteins %>% column_to_rownames("Row.names")
pd1_proteins$combo <- paste(pd1_proteins$condition, pd1_proteins$time, pd1_proteins$stimulation, sep = "_")

## CV of proteins per condition ####
proteins_cv <- aggregate(.~ combo, data = pd1_proteins[,c(1:59,66)], 
                         function(x) sd(x) / mean(x) * 100)

proteins_cv_melt <- reshape2::melt(proteins_cv)
proteins_cv_melt$proteins <- "PD1/NF-kappaB"

pdf("/home/degan/dia_mass_spec/figures/coef_variation_condition.pdf", height = 5, width = 4)
ggplot(data = proteins_cv_melt, aes(x=reorder(combo, value), y=value, color = combo)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + scale_color_npg() + labs(color = "Condition") + xlab("")
dev.off()

mean(proteins_cv_melt$value)

## random proteins ####
################################################################################
random_mean <- c()
for (i in seq(1,100,1)) {

  rand_proteins <- imputed_data[sample(nrow(imputed_data), 59),]
  rand_proteins <- t(rand_proteins)
  rand_proteins <- merge(rand_proteins, annotation_df, by =0)
  rand_proteins <- rand_proteins %>% column_to_rownames("Row.names")
  rand_proteins$combo <- paste(rand_proteins$condition, rand_proteins$time, rand_proteins$stimulation, sep = "_")

  ## CV of proteins per condition ####
  proteins_cv_rand <- aggregate(.~ combo, data = rand_proteins[,c(1:59,66)], 
                         function(x) sd(x) / mean(x) * 100)

  proteins_cv_melt_rand <- reshape2::melt(proteins_cv_rand)
  proteins_cv_melt_rand$proteins <- "Random"

  random_mean <- append(random_mean, mean(proteins_cv_melt_rand$value), i)
}

ggplot(data = proteins_cv_melt_rand, aes(x=reorder(combo, value), y=value, color = combo)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + scale_color_npg() + labs(color = "Condition_R") + xlab("")


## PATHWAY ####
################################################################################

msigdbr_df <- msigdbr(species = "human", category = "C7")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

ifg_res <- pathwaysH[["GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP"]]
ifg_res <- gene_names[gene_names$Gene.names %in% ifg_res,]
ifg_res <- imputed_data[rownames(imputed_data) %in% ifg_res$Protein.IDs,]

ifg_res <- t(ifg_res)
ifg_res <- merge(ifg_res, annotation_df, by =0)
ifg_res <- ifg_res %>% column_to_rownames("Row.names")
ifg_res$combo <- paste(ifg_res$condition, ifg_res$time, ifg_res$stimulation, sep = "_")

## CV of proteins per condition ####
proteins_cv_ifg <- aggregate(.~ combo, data = ifg_res[,c(1:19,26)], 
                              function(x) sd(x) / mean(x) * 100)

proteins_cv_melt_ifg <- reshape2::melt(proteins_cv_ifg)

mean(proteins_cv_melt_ifg$value)
proteins_cv_melt_ifg$proteins <- "PD1_VS_CTRL_Public"

ggplot(data = proteins_cv_melt_ifg, aes(x=reorder(combo, value), y=value, color = combo)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + scale_color_npg() + labs(color = "Condition_R") + xlab("")

### MERGE RESULTS ####
################################################################################

CV_values <- do.call(rbind, list(proteins_cv_melt_ifg, proteins_cv_melt, proteins_cv_melt_rand))

cv_mean <- aggregate(CV_values$value, FUN=mean, 
                     by=list(condition = CV_values$combo, proteins=CV_values$proteins))

cv_mean$condition <- factor(cv_mean$condition, 
                                      levels=c("PD1_0min_unstim", "PD1_20min_stim",  "PD1_24hr_stim", "PD1_24hr_unstim",
                                               "PD1_5min_stim","Con_0min_unstim",  "Con_20min_stim",  "Con_24hr_stim", 
                                               "Con_24hr_unstim", "Con_5min_stim"))

ggplot(cv_mean, aes(x=condition, y=x, group=proteins, color=proteins)) + geom_point() +
  geom_line() + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + scale_color_npg()

ggplot(CV_values, aes(x = value, fill = proteins)) + geom_density(alpha = 0.5) +
  theme_bw() + xlim(-5,40) + scale_fill_npg()

my_comparisons <- list(c("Random", "PD1/NF-kappaB"), c("PD1_VS_CTRL_Public", "PD1/NF-kappaB"))
ggplot(cv_mean, aes(x=proteins, y=x)) + geom_boxplot() + geom_jitter(aes(color=condition)) + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)

## ABUNDANCE OF PD1 CLUSTER AT EACH TIME POINT ####
#################################################################################

pd1_proteins_pd1 <- pd1_proteins[pd1_proteins$condition == "PD1" | pd1_proteins$time == "0min" & pd1_proteins$stimulation == "stim",]

pd1_proteins$time_condition <- paste(pd1_proteins$time, pd1_proteins$condition, sep="_")

proteins_time <- aggregate(.~ combo, data = pd1_proteins[,c(1:59,66)], mean)

proteins_time <- proteins_time %>% column_to_rownames("combo") %>% t() %>% reshape2::melt()
proteins_time$condition <- sub("_.*", "", proteins_time$Var2)

proteins_time$condition <- sub(".*_", "", proteins_time$Var2)
proteins_time$time <- factor(proteins_time$time, levels=c("0min", "5min","20min", "24hr"))


ggplot(proteins_time, aes(x=time, y=value, color = condition)) + geom_point() +
  geom_line() + theme_bw() + theme(axis.text.x = element_text(angle = 90)) 
