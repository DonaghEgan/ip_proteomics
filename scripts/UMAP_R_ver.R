## load library ####
################################################################################

library(monocle3)
library(dplyr)
library(viridis)
library(reshape2)
library(ggplot2)
library(data.table)
library(GOfuncR)
library(ggpubr)
library(stringr)

## Read in imputed data matrix and anntotations ####
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

## Store data in a cell_data_set object ####
################################################################################

cds <- new_cell_data_set(as.matrix(t(imputed_data)),
                         cell_metadata = NULL)

cds = preprocess_cds(cds, num_dim = 100, method = "PCA", norm_method = "none")

cds = reduce_dimension(cds,
                       preprocess_method = "PCA",
                       reduction_method = "UMAP", 
                       max_components = 2, 
                       umap.min_dist = 0.05, 
                       umap.n_neighbors = 10L,
                       umap.metric = "cosine")

cds = cluster_cells(cds, resolution = 1e-3, k = 5)

# plot strain umap
plot_cells(cds, 
           cell_size = 0.75,
           label_groups_by_cluster = F,
           label_cell_groups = T, 
           group_label_size = 4)

plot_cells(cds, genes=c("P16401"))

marker_test_res <- top_markers(cds, group_cells_by="cluster", cores=8)

## 13 = pdcd1 
## 1 = PTPN11

cluster_assignment <- data.frame(cluster = cds@clusters@listData[["UMAP"]][["clusters"]])
cluster_1 <- cluster_assignment[cluster_assignment$cluster == 7, , drop=F]
cluster_1 <- merge(cluster_1, pgroups.meta, by.x =0 , by.y = "Protein.IDs")

## data such as functional and notation and gene ontology ######################

## create input dataframe with candidate and background genes
candi_gene_ids = cluster_1$Gene.names
candi_gene_ids <- strsplit(candi_gene_ids, ";", fixed=F)
candi_gene_ids <- Reduce(c,candi_gene_ids)

bg_gene_ids = pgroups.meta$Gene.names[!pgroups.meta$Gene.names %in% candi_gene_ids]
bg_gene_ids <- strsplit(bg_gene_ids, ";", fixed=F)
bg_gene_ids <- Reduce(c,bg_gene_ids)

is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),
                            is_candidate)

res_hyper_bg = go_enrich(input_hyper_bg, n_randsets=1000)

top_gos_hyper = res_hyper_bg[[1]][1:5, 'node_id']
plot_anno_scores(res_hyper_bg, top_gos_hyper)

##
###############################################################################

## cluster 13 pd1 
genes_cluster_pd1 <- cluster_1$Row.names
genes_cluster_pd1 <- imputed_data[which(rownames(imputed_data) %in% genes_cluster_pd1), colnames(imputed_data) %in% rownames(annotation_df)]
genes_cluster_pd1 <- t(genes_cluster_pd1)
annotation_df <- annotation_df[annotation_df$condition == "PD1", ]
genes_cluster_pd1_annot <- cbind(genes_cluster_pd1, annotation_df)
avg_cluster_pd1 <- df <- aggregate(.~ antibody, genes_cluster_pd1_annot[,c(1:59,61)], mean)
avg_cluster_pd1 <- avg_cluster_pd1 %>% column_to_rownames("antibody")
avg_cluster_pd1 <- t(avg_cluster_pd1)
avg_cluster_pd1 <- avg_cluster_pd1[which(rownames(avg_cluster_pd1) != "Q15116"),]
avg_cluster_pd1 <- melt(avg_cluster_pd1)
## remove pd1
comparisons <- list( c("Csa", "Niv"), c("Csa", "Pem"), c("Niv", "Pem") )
ggplot(avg_cluster_pd1, aes(x=Var2, y=value)) + 
  geom_boxplot() +  stat_compare_means(method = "wilcox.test", comparisons = comparisons) + theme_classic()


## cluster 7 TCR 
genes_cluster_pd1 <- cluster_1$Row.names
genes_cluster_pd1 <- imputed_data[which(rownames(imputed_data) %in% genes_cluster_pd1),]
genes_cluster_pd1 <- t(genes_cluster_pd1)
genes_cluster_pd1_annot <- cbind(genes_cluster_pd1, annotation_df)
avg_cluster_pd1 <- df <- aggregate(.~ condition, genes_cluster_pd1_annot[,c(1:118)], mean)
avg_cluster_pd1 <- avg_cluster_pd1 %>% column_to_rownames("condition")
avg_cluster_pd1 <- t(avg_cluster_pd1)
avg_cluster_pd1 <- avg_cluster_pd1[which(rownames(avg_cluster_pd1) != "Q15116"),]
avg_cluster_pd1 <- melt(avg_cluster_pd1)
## remove pd1
ggplot(avg_cluster_pd1, aes(x=Var2, y=value)) + 
  geom_boxplot() +  stat_compare_means() + theme_bw()


## turn clusters to a list to be scored ####
cluster_assignment <- data.frame(cluster = cds@clusters@listData[["UMAP"]][["clusters"]])
cluster_assignment$protein <- rownames(cluster_assignment)
cluster_assignment <- cluster_assignment  %>% group_by(cluster) %>% summarise(across(everything(), str_c, collapse=" ")) 
cluster_assignment$protein <- strsplit(cluster_assignment$protein, split = " ")
cluster_list <- cluster_assignment$protein

## score list ssgsea ####
library(GSVA)

clusters_scored <-  gsva(as.matrix(imputed_data), cluster_list, mx.diff=FALSE, verbose=TRUE, parallel.sz=1, method = "ssgsea") #ssgsea
rown
pheatmap(clusters_scored, show_colnames = F, scale = "row", annotation_col = annotation_df, show_rownames = T)


####
diff_exp_protein <- diff_exp_protein[diff_exp_protein$type == "PD1vsCtr" & diff_exp_protein$logFC > 0 &
                                       diff_exp_protein$adj.P.Val < 0.05,]

intersect(diff_exp_protein$rn, cluster_1$Row.names)



