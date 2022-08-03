diff_exp_filter <- readRDS("/home/degan/ip_proteomics/inputs/diff_bound_filter.Rds")
cluster_info <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
gene_names <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

for (i in unique(cluster_info$cluster_id)) {
  genes_cluster <- cluster_info[cluster_info$cluster_id == i, , drop=F]
  genes_cluster <- imputed_data[which(rownames(imputed_data) %in% genes_cluster$Row.names), ]
  genes_cluster_annot <- cbind(t(genes_cluster), annotation_df)
  
  avg_cluster <- aggregate(.~ condition, genes_cluster_annot[,c(1:(ncol(genes_cluster_annot) - 5))], mean)
  avg_cluster <- avg_cluster %>% column_to_rownames("condition")
  avg_cluster <- reshape2::melt(t(avg_cluster))
  
  temp_name <- i
  temp_name <- gsub("/", "_", temp_name)
  
  ggplot(avg_cluster, aes(x=Var2, y=value, fill=Var2)) + 
  geom_boxplot(alpha = 0.7, width = 0.3)  + geom_jitter(width = 0.1, alpha = 0.5) +
  stat_compare_means(method = "t.test") + theme_classic() + xlab("Condition") +
  ylab("Avg protein Abundance") +  theme(legend.title=element_blank(),
                                         legend.position="none") + scale_fill_npg() + ggtitle(temp_name)
  ggsave(paste(paste("/home/degan/ip_proteomics/figures/UMAP/avg_exp_pd1_con/", temp_name, sep = "cluster_"), ".pdf", sep=""), height = 3, width = 2.5)

}

## clusters statistically diff:
## TCR complex, Mitochondrial ribosome, Nuclear Lumen, Mito large ribo sub, CCR4−NOT core complex,
## PD1_NF−kappaB, Nuclear Protein comp, MLL 1_2, CCR4−NOT complex, Ribonucleoprotein granule, 
## MCM−CMG complex


