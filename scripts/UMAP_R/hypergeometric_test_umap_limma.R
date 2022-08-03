## Performing hypergeometric test between complexes and DAP ####
################################################################################

library(corrplot)
library("gplots")

diff_exp_filter <- readRDS("/home/degan/ip_proteomics/inputs/diff_bound_filter.Rds")
cluster_info <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
gene_names <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds"))


## Turn differentially abundant proteins between conditions into a list ###
gene_list <- list("PD1.0min" = filter(diff_exp_filter, (type == "PD1_0vsCon_0" & logFC > 0)),
                  "PD1.5min" = filter(diff_exp_filter, (type == "PD1_5vsCon_5" & logFC > 0)),
                  "PD1.20min" = filter(diff_exp_filter, (type == "PD1_20vsCon_20" & logFC > 0)),
                  "PD1.24hr" = filter(diff_exp_filter, (type == "PD1_24vsCon_24" & logFC > 0)),
                  "PD1.time" = filter(diff_exp_filter, (type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20"))),
                  "Ctr.0min" = filter(diff_exp_filter, (type == "PD1_0vsCon_0" & logFC < 0)),
                  "Ctr.5min" = filter(diff_exp_filter, (type == "PD1_5vsCon_5" & logFC < 0)),
                  "Ctr.20min" = filter(diff_exp_filter, (type == "PD1_20vsCon_20" & logFC < 0)),
                  "Ctr.24hr" = filter(diff_exp_filter, (type == "PD1_24vsCon_24" & logFC < 0)))

## Turn clusters to a list ####
cluster_ids <- cluster_info  %>% group_by(cluster_id) %>% summarise(across(everything(), str_c, collapse=" ")) 
cluster_ids$Row.names <- strsplit(cluster_ids$Row.names, split = " ")
cluster_list <- cluster_ids$Row.names
names(cluster_list) <- cluster_ids$cluster_id

## test overlap ####
PopSize <- nrow(imputed_data)
hyper_test <- list()
for (complex in names(cluster_list)) {
  hyper_complex <- list()
  for (condition in names(gene_list)) {
    overlap = length(intersect(cluster_list[[complex]], gene_list[[condition]][["rn"]]))
    list1 <- length(gene_list[[condition]][["rn"]])
    list2 <- length(cluster_list[[complex]])
    pval = phyper(overlap-1, list1, PopSize-list1, list2, lower.tail = FALSE, log.p = FALSE)
    hyper_complex[[condition]] <- pval
  }
  hyper_test[[complex]] <- hyper_complex
  padjust <- p.adjust(hyper_test[[complex]], method = "BH")
  hyper_test[[complex]] <- -log(padjust)
}

## visualize results ####
hyper_test_df <- data.frame(matrix(unlist(hyper_test), nrow=length(hyper_test), byrow=TRUE),
                            row.names = names(hyper_test))
colnames(hyper_test_df) <- names(hyper_test[["MLL 1/2"]])
pheatmap(hyper_test_df, cluster_rows = T, cluster_cols = T)
corrplot(as.matrix(hyper_test_df), is.corr = FALSE, tl.cex = .7)

hyper_test_df_m <- reshape2::melt(as.matrix(hyper_test_df))
hyper_test_df_m$adjpval <- as.factor(ifelse(hyper_test_df_m$value > 1.3010, "Significant", "Not-significant")) 

pdf("/home/degan/ip_proteomics/figures/UMAP/DBP_complexes_ass.pdf", width = 5, height = 4)
ggplot(hyper_test_df_m, aes(x= Var2, y=Var1, size=value, color=adjpval)) + 
  geom_point(alpha = 0.8) + scale_color_npg() +
  theme_classic() + theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  labs(size='-log(padjust)') + ylab("Protein Cluster") + xlab("Condition")
dev.off()


tcr_pd1_0 <- data.frame(inter = intersect(cluster_list[["TCR complex"]], gene_list[["PD1.0min"]][["rn"]]))
tcr_pd1_0 <- merge(tcr_pd1_0, gene_names, by.x="inter", by.y="Protein.IDs")

tcr_cntr_0 <- data.frame(inter = intersect(cluster_list[["TCR complex"]], gene_list[["Ctr.0min"]][["rn"]]))
tcr_cntr_0 <- merge(tcr_cntr_0, gene_names, by.x="inter", by.y="Protein.IDs")


