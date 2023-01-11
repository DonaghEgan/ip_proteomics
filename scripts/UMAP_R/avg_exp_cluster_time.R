library(ggsci)

cluster_info <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
gene_names <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

## PD1 samples only
annotation_pd1 <- annotation_df[annotation_df$type == "PD1",]
imputed_pd1 <- imputed_data[,colnames(imputed_data) %in% rownames(annotation_pd1)]

plot_list <- list()
for (i in unique(cluster_info$cluster)) {
  genes_cluster <- cluster_info[cluster_info$cluster == i, , drop=F]
  genes_cluster <- imputed_data[which(rownames(imputed_pd1) %in% genes_cluster$Row.names), ]
  genes_cluster_annot <- cbind(t(genes_cluster), annotation_pd1)
  
  avg_cluster <- aggregate(.~ time, genes_cluster_annot[ ,-which(names(genes_cluster_annot) %in%
                                                                   c("type","antibody", "batch","stimulation","condition"))], mean)
  avg_cluster <- avg_cluster %>% column_to_rownames("time")
  avg_cluster <- reshape2::melt(t(avg_cluster))
  avg_cluster$Var2 <- factor(avg_cluster$Var2, levels=c("0min", "5min", "20min", "24hr"))
  
  temp_name <- i
  temp_name <- gsub("/", "_", temp_name)
  
  p1 <- ggplot(avg_cluster, aes(x=Var2, y=value, fill=Var2)) + 
    geom_boxplot(alpha = 0.7, width = 0.3)  + geom_jitter(width = 0.1, alpha = 0.5, size=0.1) +
    theme_classic() + theme(legend.title=element_blank(),
                            legend.position="none", plot.title=element_text(margin=margin(t=10,b=-5), size=9)) + scale_fill_npg() + 
    labs(
      title = temp_name,
      y = "", x = ""
    ) +
    stat_compare_means(method = "anova", size=2.5)
  plot_list[[i]] <- p1
}

## PTPN11 cluster per antibody 
cluster_19 <- cluster_info[cluster_info$cluster ==19,]
cluster_19 <- imputed_data[rownames(imputed_data) %in% cluster_19$Row.names,]
cluster_19 <- data.frame(t(cluster_19))
cluster_19 <- merge(cluster_19, annotation_df[,c(6), drop=F], by=0)
cluster_19 <- column_to_rownames(cluster_19, "Row.names")
cluster_19 <- aggregate(.~ condition, cluster_19, mean)
cluster_19$condition <- substr(cluster_19$condition, 5,7)
cluster_19_melt <- reshape2::melt(cluster_19, "condition")

my_comparisons <- list( c("Csa", "Niv"), c("Csa", "Pem"), c("Niv", "Pem"))
pdf("/home/degan/ip_proteomics/figures/UMAP/cluster_shp2_antibody.pdf", height = 4, width = 4)
ggplot(cluster_19_melt, aes(x=condition, y=value, fill="#e64b35")) + geom_boxplot() + geom_point() + 
  stat_compare_means(method="wilcox.test", comparisons = my_comparisons) + theme_bw() + theme(legend.position = "None")
dev.off()
  

