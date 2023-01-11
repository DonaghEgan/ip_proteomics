################################################################################


## load library ####
################################################################################

library(monocle3)
library(ggsci)
library(dplyr)
library(viridis)
library(reshape2)
library(ggplot2)
library(data.table)
library(GOfuncR)
library(ggpubr)
library(stringr)
library(lme4)
library(readxl)
library(gdata)
library(monocle3)
library(gprofiler2)
library(RColorBrewer)
library(cluster)

## Read in imputed data matrix and anntotations ####
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")
pgroups.meta <- readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds")

## Store data in a cell_data_set object ####
################################################################################

set.seed(16351893)

cds <- new_cell_data_set(as.matrix(t(imputed_data)),
                         cell_metadata = NULL)

cds = preprocess_cds(cds, num_dim = 80, method = "PCA", norm_method = "none", scaling = T)

cds = reduce_dimension(cds,
                       preprocess_method = "PCA",
                       reduction_method = "UMAP", 
                       max_components = 2, 
                       umap.min_dist = 0.05, 
                       umap.n_neighbors = 10L,
                       umap.metric = "cosine")

## cluster cells ####
################################################################################
sil_stats <- list()
for (i in seq(0.000001, 0.005,  0.0001)) {
  print(i)
  cds = cluster_cells(cds, resolution = i, k = 5)
  
  cluster_assignment <- data.frame(cluster = cds@clusters@listData[["UMAP"]][["clusters"]])
  cluster_assignment <- merge(cluster_assignment, pgroups.meta, by.x =0 , by.y = "Protein.IDs")
  
  umap_embedding <- cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]]
  umap_embedding <- merge(umap_embedding, cluster_assignment, by.x=0, by.y="Row.names")
  dist.matrix <- dist(x = umap_embedding[,c(2:3)])

  sil <- silhouette(x = as.numeric(x = as.factor(x = umap_embedding$cluster)), dist = dist.matrix)
  sil_sum <- summary(sil)
  mean_sil_sum <- mean(sil_sum$avg.width)
  sil_stats[[as.character(i)]] <- c(mean_sil_sum)
}

cds = cluster_cells(cds, resolution = 0.000501, k = 5)

plot_cells(cds, 
           cell_size = 0.75,
           label_groups_by_cluster = F,
           label_cell_groups = T, 
           group_label_size = 4)

cluster_assignment <- data.frame(cluster = cds@clusters@listData[["UMAP"]][["clusters"]])
cluster_assignment <- merge(cluster_assignment, pgroups.meta, by.x =0 , by.y = "Protein.IDs")


sil.df <- do.call(rbind.data.frame, sil_stats)
sil.df$resolution <- as.numeric(names(sil_stats))
colnames(sil.df)[1] <- "Mean SI score"

pdf("/home/degan/ip_proteomics/figures/UMAP/SI_index_clustering.pdf", height = 4, width = 4)
ggplot(sil.df, aes(x=resolution, y=`Mean SI score`)) + geom_point() + 
  geom_line(color="#5fbcd3", alpha=0.5) + theme_bw() + xlab("Resolution Parameter")
dev.off()
