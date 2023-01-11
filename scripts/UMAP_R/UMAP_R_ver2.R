################################################################################
# Donagh Egan
# POIAZ 
# Date: 30nd June 2022
# 

# Description: Performs UMAP clustering on the ip proteomics; labels the resulting 
# clusters using GO; and analyzes the associations between clusters and outcomes.
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

## Read in imputed data matrix and anntotations ####
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")
pgroups.meta <- readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds")

## Store data in a cell_data_set object ####
################################################################################

set.seed(16351893)

cds <- new_cell_data_set(as.matrix(t(imputed_data)),
                         cell_metadata = NULL)

cds = preprocess_cds(cds, num_dim = 80, method = "PCA", norm_method = "none")

cds = reduce_dimension(cds,
                       preprocess_method = "PCA",
                       reduction_method = "UMAP", 
                       max_components = 2, 
                       umap.min_dist = 0.05, 
                       umap.n_neighbors = 10L,
                       umap.metric = "cosine")

## cluster cells ####
################################################################################

cds = cluster_cells(cds, resolution = 0.000501, k = 5)

## plot UMAP ####
################################################################################

pdf("/home/degan/ip_proteomics/figures/UMAP/UMAP_clusters.pdf", height = 4, width = 4)
plot_cells(cds, 
           cell_size = 0.5,
           label_groups_by_cluster = F,
           label_cell_groups = T, 
           group_label_size = 4)
dev.off()

## protein cluster assignment ####
################################################################################

cluster_assignment <- data.frame(cluster = cds@clusters@listData[["UMAP"]][["clusters"]])
cluster_assignment <- merge(cluster_assignment, pgroups.meta, by.x =0 , by.y = "Protein.IDs")

## GO enrichment - each cluster ####
################################################################################
GO_results <- list()
for (i in unique(cluster_assignment$cluster)) {
  cluster_proteins <- cluster_assignment[cluster_assignment$cluster == i, , drop=F]
  
  ## create input dataframe with candidate and background genes
  candi_gene_ids = cluster_proteins$Gene.names
  candi_gene_ids <- strsplit(candi_gene_ids, ";", fixed=F)
  candi_gene_ids <- Reduce(c,candi_gene_ids)
  
  bg_gene_ids = pgroups.meta$Gene.names[!pgroups.meta$Gene.names %in% candi_gene_ids]
  bg_gene_ids <- strsplit(bg_gene_ids, ";", fixed=F)
  bg_gene_ids <- Reduce(c,bg_gene_ids)
  
  is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
  input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),
                              is_candidate)
  
  res_hyper_bg = go_enrich(input_hyper_bg, n_randsets=1000)
  
  ## store in list 
  GO_results[[i]] <- res_hyper_bg
}

## profiler - includes kegg etc ####
################################################################################

gprofiler_results <- list()
for (i in unique(cluster_assignment$cluster)) {
  cluster_proteins <- cluster_assignment[cluster_assignment$cluster == i, , drop=F]
  
  ## create input dataframe with candidate and background genes
  candi_gene_ids = cluster_proteins$Gene.names
  candi_gene_ids <- strsplit(candi_gene_ids, ";", fixed=F)
  candi_gene_ids <- Reduce(c,candi_gene_ids)
  
  gostres <- gost(query = candi_gene_ids,correction_method = c("fdr"), 
                  organism = "hsapiens", significant = T, user_threshold = 0.05, 
                  sources = c("GO", "KEGG", "REAC"))
  
  ## store in list 
  gprofiler_results[[i]] <- gostres
}


## Annotate Clusters according to GO ####
################################################################################

colData(cds)$cluster_id <- cds@clusters@listData[["UMAP"]][["clusters"]]

colData(cds)$cluster_id <- dplyr::recode(colData(cds)$cluster_id,
                                                 "1"="MLL 1/2",
                                                 "2"="U2 spliceosomal complex",
                                                 "3"="Mitochondrial ribosome",
                                                 "4"="SWI/SNF complex",
                                                 "5"="Nuclear Lumen",
                                                 "6"="ER",
                                                 "7"="TCR complex",
                                                 "8"="Large Ribosomal subunit",
                                                 "9"="Cytosolic ribosome",
                                                 "10"="CCR4-NOT complex",
                                                 "11"="CCR4-NOT core complex",
                                                 "12"="Ribonucleoprotein granule",
                                                 "13"="PD1/NF-kappaB",
                                                 "14"="Mito large ribo sub",
                                                 "15"="48/43S pre-initiation comp",
                                                 "16"="Nuclear Protein comp",
                                                 "17"="MCM-CMG complex",
                                                 "18"="DNA repair",
                                                 "19"="PI3K",
                                                 "20"="PIP2")

plot_cells(cds, group_cells_by="cluster", color_cells_by="cluster_id",
           cell_size = 0.75, group_label_size = 3)

## Avg Exp of PD1 and TCR in ctrl and pd1 ####
################################################################################

## PD1
genes_cluster_pd1 <- cluster_assignment[cluster_assignment$cluster == 13, , drop=F]
genes_cluster_pd1 <- imputed_data[which(rownames(imputed_data) %in% genes_cluster_pd1$Row.names), ]
genes_cluster_pd1_annot <- cbind(t(genes_cluster_pd1), annotation_df)

## average of each gene per condition (PD1/Ctrl)
avg_cluster_pd1 <- aggregate(.~ type, genes_cluster_pd1_annot[,c(1:69)], mean)
avg_cluster_pd1 <- avg_cluster_pd1 %>% column_to_rownames("type")
avg_cluster_pd1 <- reshape2::melt(t(avg_cluster_pd1))

## draw and save plot #### 
pdf("/home/degan/ip_proteomics/figures/UMAP/pd1_cluster_conVSpd1.pdf", height = 4, width = 4)
p1 <- ggplot(avg_cluster_pd1, aes(x=Var2, y=value, fill=Var2)) + 
  geom_boxplot() +  stat_compare_means(method = "t.test") + theme_classic() + xlab("Condition") +
  ylab("Mean PD1/NF-kappaB Abd") +  theme(legend.title=element_blank()) + scale_fill_npg()
dev.off()

## average of each gene per antibody (NIV/PEM/CSA)
avg_cluster_ant <- aggregate(.~ antibody, genes_cluster_pd1_annot[,c(1:59, 61)], mean)
avg_cluster_ant <- avg_cluster_ant %>% column_to_rownames("antibody")
avg_cluster_ant<- melt(t(avg_cluster_ant))

## draw and save plot #### 
my_comparisons <- list( c("Niv", "Csa"), c("Niv", "Pem"), c("Pem", "Csa"))
pdf("/home/degan/ip_proteomics/figures/UMAP/pd_cluster_ant_labelled.pdf", height = 4, width = 4)
ggplot(avg_cluster_ant, aes(x=Var2, y=value, fill=Var2, label=Var1)) + 
  geom_boxplot(alpha = 0.7) + geom_point(alpha = 0.5) +  
  stat_compare_means(method = "t.test", comparisons = my_comparisons) + theme_classic() + xlab("Antibody") +
  ylab("PD1/NF-kappaB Abundance") +  theme(legend.title=element_blank()) + scale_fill_npg() + geom_text(aes(label=ifelse(Var1 == "Q15116",as.character(Var1),'')), 
                                                                                                        hjust=0, vjust=0)
dev.off()

## TCR
genes_cluster_tcr <- cluster_assignment[cluster_assignment$cluster == 17, , drop=F]
genes_cluster_tcr <- imputed_data[which(rownames(imputed_data) %in% genes_cluster_tcr$Row.names), ]
genes_cluster_tcr_annot <- cbind(t(genes_cluster_tcr), annotation_df)

## average of each gene per antibody (csa,niv,pem)
avg_cluster_tcr <- df <- aggregate(.~ type, genes_cluster_tcr_annot[,c(1:10)], mean)
avg_cluster_tcr <- avg_cluster_tcr %>% column_to_rownames("type")
avg_cluster_tcr <- melt(t(avg_cluster_tcr))

## draw and save plot 
pdf("/home/degan/ip_proteomics/figures/UMAP/tcr_cluster_conVSpd1.pdf", height = 4, width = 4)
p2 <- ggplot(avg_cluster_tcr, aes(x=Var2, y=value, fill=Var2)) + 
  geom_boxplot(alpha = 0.7)  +
  stat_compare_means(method = "t.test") + theme_classic() + xlab("Condition") +
  ylab("TCR Complex") +  theme(legend.title=element_blank())
dev.off()

## average of each gene per antibody (NIV/PEM/CSA)
avg_cluster_ant <- df <- aggregate(.~ antibody, genes_cluster_pd1_annot[,c(1:59, 61)], mean)
avg_cluster_ant <- avg_cluster_ant %>% column_to_rownames("antibody")
avg_cluster_ant<- melt(t(avg_cluster_ant))

## Check Correlation between complexes: cluster ####
################################################################################

## turn clusters to a list to be scored ####
cluster_ids <- data.frame(cluster_id =  cds@clusters@listData[["UMAP"]][["clusters"]])
cluster_ids$protein <- rownames(cluster_ids)
cluster_ids <- cluster_ids  %>% group_by(cluster_id) %>% summarise(across(everything(), str_c, collapse=" ")) 
cluster_ids$protein <- strsplit(cluster_ids$protein, split = " ")
cluster_list <- cluster_ids$protein
names(cluster_list) <- cluster_ids$cluster_id

## score list ssgsea ####
library(GSVA)
clusters_scored <-  gsva(as.matrix(imputed_data), cluster_list, mx.diff=FALSE, verbose=TRUE, parallel.sz=1, method = "ssgsea") #ssgsea
clusters_cor <- cor(t(clusters_scored), method = "spearman")
pdf("/home/degan/ip_proteomics/figures/UMAP/cluster_correlations.pdf", height = 4, width = 5)
pheatmap(clusters_cor, show_colnames = F, fontsize_row = 8, treeheight_row = 14, 
         treeheight_col = 14, show_rownames = T)
dev.off()

## Correlation tcr and pd1 ####
################################################################################

scores_pd1_tcr <- data.frame(t(clusters_scored[c("PD1/NF-kappaB","TCR complex"),]))
scores_pd1_tcr <- cbind(scores_pd1_tcr, annotation_df)

## plotting relationship TCR and PD1 according to stim ####
pdf("/home/degan/ip_proteomics/figures/UMAP/tcr_pd1_plot_fill_stim.pdf", height = 4, width = 4)
ggplot(scores_pd1_tcr , aes(x=PD1.NF.kappaB, y=TCR.complex)) + geom_point(aes(color=stimulation)) +
       geom_smooth(method=lm, se=FALSE, color = "black",
                   linetype="dashed") + theme_classic()
dev.off()

## relationship mcm and pd1 ####

scores_pd1_mcm <- data.frame(t(clusters_scored[c("PD1/NF-kappaB","MCM-CMG complex"),]))
scores_pd1_mcm <- cbind(scores_pd1_mcm, annotation_df)

pdf("/home/degan/ip_proteomics/figures/UMAP/PD1_MCM_antibody.pdf", height = 5, width = 5)
ggplot(scores_pd1_mcm , aes(x=PD1.NF.kappaB, y=MCM.CMG.complex, color = antibody)) + geom_point() +
  geom_smooth(method=lm, se=TRUE, linetype="dashed") + theme_classic()
dev.off()

source("scripts/t_test_reg.R", chdir = TRUE)
model_csa <- scores_pd1_mcm[which(scores_pd1_mcm$antibody == "Csa"),]
model_niv <- scores_pd1_mcm[which(scores_pd1_mcm$antibody == "Niv"),]
model_pem <- scores_pd1_mcm[which(scores_pd1_mcm$antibody == "Pem"),]

res_niv_csa <- ttest_reg(y1 = model_csa$MCM.CMG.complex, x1 = model_csa$PD1.NF.kappaB,
                         y2 = model_niv$MCM.CMG.complex, x2 = model_niv$PD1.NF.kappaB,
                         y3 = model_pem$MCM.CMG.complex, x3 = model_pem$PD1.NF.kappaB)

res_niv_pem <- ttest_reg(y1 = model_pem$PD1.NF.kappaB, x1 = model_pem$MCM.CMG.complex,
                         y2 = model_niv$PD1.NF.kappaB, x2 = model_niv$MCM.CMG.complex)

res_csa_pem <- ttest_reg(y1 = model_pem$PD1.NF.kappaB, x1 = model_pem$MCM.CMG.complex,
                         y2 = model_csa$PD1.NF.kappaB, x2 = model_csa$MCM.CMG.complex)


test <- ttest_reg(y1 = model_pem$PD1.NF.kappaB, x1 = model_pem$MCM.CMG.complex,
                  y2 = model_pem$PD1.NF.kappaB, x2 = model_pem$MCM.CMG.complex)

reg1 = glm(model_csa$PD1.NF.kappaB~model_csa$MCM.CMG.complex)
reg2 = glm(model_niv$PD1.NF.kappaB~model_csa$MCM.CMG.complex)

#ant_mixed = lmer(MCM.CMG.complex ~ PD1.NF.kappaB + (1 +  PD1.NF.kappaB| antibody), data = scores_pd1_mcm)

## Overlay with CRAPOME
################################################################################

CrapDB <- read_csv("/home/degan/ip_proteomics/inputs/crapome/crap_db.csv")
CrapDB <- CrapDB[CrapDB$NUM_EXPT >= 5,]
crap_proteins <- CrapDB$GENE

nf_kb_go <- read.xlsx("/home/degan/ip_proteomics/inputs/GO/GO_term_summary_20220630_090345.xlsx")
nf_kb_go <- toupper(nf_kb_go$Symbol)

gene_names_ip <- pgroups.meta$Gene.names
gene_names_ip <- strsplit(gene_names_ip, ";", fixed=F)
gene_names_ip <- Reduce(c,gene_names_ip)

crap_ip <- intersect(crap_proteins, gene_names_ip)
is_crap <- c(rep(1,length(crap_ip)))
not_crap <- c(rep(0, length(gene_names_ip[!gene_names_ip %in% crap_ip])))

crap_df = data.frame(Gene.names = c(crap_ip, gene_names_ip[!gene_names_ip %in% crap_ip]),
                           crap = c(is_crap, not_crap))

cluster_assignment <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
cluster_assignment <- cluster_assignment[cluster_assignment$cluster %in% c(8,9),]
cluster_crap <- merge(crap_df, cluster_assignment)
saveRDS(cluster_crap, "/mnt/Data/Vadim/POIAZ/Donagh/proteomics22_unified/downstream_analysis/saved_objects/ribo_cluster_CRAPOME.Rds")

crap_df <- merge(crap_df, pgroups.meta, by="Gene.names")
crap_df <- crap_df[crap_df$crap == 1,]

colData(cds)$crap_proteins <- cds@clusters@listData[["UMAP"]][["clusters"]]
for (i in names(colData(cds)$crap_proteins)) {
    if (i %in% crap_df$Protein.IDs) {
    colData(cds)$crap_proteins[[i]] <- 1 
   } else {
    colData(cds)$crap_proteins[[i]] <- 2
   }
}

pdf("/home/degan/ip_proteomics/figures/UMAP/umap_crapome_annotation.pdf", height = 5, width = 5)
plot_cells(cds, group_cells_by="cluster", color_cells_by="crap_proteins",
           cell_size = 1, group_label_size = 3)  
dev.off()

## plot GO results ####
#################################################################################

cluster_pd1 <- data.frame(GO_results[["13"]])
nf_terms <- c("NIK/NF-kappaB signaling", "regulation of NIK/NF-kappaB signaling")
cluster_pd1 <- rbind(cluster_pd1[cluster_pd1$node_name %in% nf_terms, ],
                     cluster_pd1[c(69:72), ])
cluster_pd1$logpval <- -log(cluster_pd1$FWER_overrep)

pdf("/home/degan/ip_proteomics/figures/UMAP/GO_term_pd1.pdf", height = 4, width = 5)
ggplot(data=cluster_pd1, aes(x=node_name, y=logpval)) +
  geom_bar(stat="identity", fill = "blue") + coord_flip() + xlab("GO Term") + ylab("-Log(adj Pal)") + 
  theme_bw()
dev.off()

## TCR ####
cluster_tcr <- data.frame(GO_results[["7"]])
cluster_tcr  <- cluster_tcr[c(1:5),]
cluster_tcr$logpval <- -log(cluster_tcr$raw_p_overrep)

pdf("/home/degan/ip_proteomics/figures/UMAP/GO_term_tcr.pdf", height = 4, width = 5)
ggplot(data=cluster_tcr, aes(x=node_name, y=logpval)) +
  geom_bar(stat="identity", fill = "red") + coord_flip() + xlab("GO Term") + ylab("-Log(raw Pval)") + 
  theme_bw()
dev.off()

## save cluster assignment df ####
saveRDS(cluster_assignment, "/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")

## prop pd1 proteins up in ctrl vs pd1 ####
pd1_proteins <- cluster_assignment[cluster_assignment$cluster == 13,]


## checking prop of pd1 cluster in dep of pd1 vs ctrl ####
################################################################################

diff_exp_protein <- diff_exp_protein[diff_exp_protein$type == "PD1vsCtr",]

pd1_up <- diff_exp_protein[diff_exp_protein$logFC > .5 & diff_exp_protein$adj.P.Val < 0.05, ]
ctrl_up <- diff_exp_protein[diff_exp_protein$logFC < -.5 & diff_exp_protein$adj.P.Val < 0.05]

prop_pd1_up <- (length(intersect(pd1_proteins$Row.names, pd1_up$rn)) / length(pd1_up$rn)) * 100
prop_ctrl_up <- (length(intersect(pd1_proteins$Row.names, ctrl_up$rn)) / length(ctrl_up$rn)) * 100


pd_up_cluster <- intersect(pd1_proteins$Row.names, pd1_up$rn)

library(ggplot2)

pd1_prop_df <- melt(data.frame(proteins_pd1_cluster = length(intersect(pd1_proteins$Row.names, pd1_up$rn)),
                          upregulated_pd1 = length(pd1_up$rn) - length(intersect(pd1_proteins$Row.names, pd1_up$rn))))

ctrl_prop_df <- melt(data.frame(proteins_pd1_cluster = length(intersect(pd1_proteins$Row.names, ctrl_up$rn)),
                               upregulated_ctrl = length(ctrl_up$rn) - length(intersect(pd1_proteins$Row.names, ctrl_up$rn))))

# Basic piechart
pdf("/home/degan/ip_proteomics/figures/UMAP/prop_pd1cluster_uppd1.pdf", width = 4, height = 4)
ggplot(pd1_prop_df, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)  + theme_void()
dev.off()
  
pdf("/home/degan/ip_proteomics/figures/UMAP/prop_pd1cluster_upctrl.pdf", width = 4, height = 4)
ggplot(ctrl_prop_df, aes(x="", y=value, fill=variable)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)  + theme_void()
dev.off()

## checking prop of pd1 cluster in dep of pd1 across time ####
################################################################################

diff_exp_protein_pd1 <- diff_exp_protein[diff_exp_protein$type %in% c("ST24vST20", "ST20vST5",
                                                                      "ST5vNsT0 "),]

diff_exp_protein_pd1 <- diff_exp_protein_pd1[abs(diff_exp_protein_pd1$logFC) > 1 & diff_exp_protein_pd1$adj.P.Val < 0.05, ]


prop_pd1_within <- (length(intersect(pd1_proteins$Row.names, diff_exp_protein_pd1$rn)) / length(diff_exp_protein_pd1$rn)) * 100



## remove crap proteins ####
################################################################################

cds <- new_cell_data_set(as.matrix(t(imputed_data[!rownames(imputed_data) %in% crap_df$Protein.IDs,])),
                         cell_metadata = NULL)

cds = preprocess_cds(cds, num_dim = 100, method = "PCA", norm_method = "none")

cds = reduce_dimension(cds,
                       preprocess_method = "PCA",
                       reduction_method = "UMAP", 
                       max_components = 2, 
                       umap.min_dist = 0.05, 
                       umap.n_neighbors = 10L,
                       umap.metric = "cosine")

cds = cluster_cells(cds, resolution = 1e-3, k = 6)

plot_cells(cds, 
           cell_size = 0.75,
           label_groups_by_cluster = F,
           label_cell_groups = T, 
           group_label_size = 4)

################################################################################

cluster_assignment <- data.frame(cluster = cds@clusters@listData[["UMAP"]][["clusters"]], 
                                 cluster_id = colData(cds)$cluster_id)
cluster_assignment <- merge(cluster_assignment, pgroups.meta, by.x =0 , by.y = "Protein.IDs")
saveRDS(cluster_assignment, "/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
tst <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")

GO_results_filter <- list()
for (i in unique(cluster_assignment$cluster)) {
  cluster_proteins <- cluster_assignment[cluster_assignment$cluster == i, , drop=F]
  
  ## create input dataframe with candidate and background genes
  candi_gene_ids = cluster_proteins$Gene.names
  candi_gene_ids <- strsplit(candi_gene_ids, ";", fixed=F)
  candi_gene_ids <- Reduce(c,candi_gene_ids)
  
  bg_gene_ids = pgroups.meta$Gene.names[!pgroups.meta$Gene.names %in% candi_gene_ids]
  bg_gene_ids <- strsplit(bg_gene_ids, ";", fixed=F)
  bg_gene_ids <- Reduce(c,bg_gene_ids)
  
  is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
  input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),
                              is_candidate)
  
  res_hyper_bg = go_enrich(input_hyper_bg, n_randsets=1000)
  
  ## subset: FDR < 0.05
  res_hyper_bg = res_hyper_bg$results
  
  ## store in list 
  GO_results_filter[[i]] <- res_hyper_bg
}


## Figure of GO term scores ####
GO_scores <- read_xlsx("/home/degan/ip_proteomics/inputs/GO/GO_IP_proteomics.xlsx")
GO_scores$logFWER <- -log(GO_scores$FWER_overrep + 0.01)

pdf("/home/degan/ip_proteomics/figures/UMAP/GO_term_each_cluster.pdf", height = 4, width = 5)
ggplot(data=GO_scores, aes(x=node_name, y=logFWER)) +
  geom_bar(stat="identity", fill = "red") + coord_flip() + xlab("GO Term") + ylab("-Log(FWER overrep)") + 
  theme_bw()
dev.off()

## Build and export graph #####

umap_graph <- cds@clusters@listData[["UMAP"]][["cluster_result"]][["relations"]]
umap_graph$from_g <- pgroups.meta$Gene.names[match(umap_graph$from,pgroups.meta$Protein.IDs)]
umap_graph$to_g <- pgroups.meta$Gene.names[match(umap_graph$to,pgroups.meta$Protein.IDs)] 
umap_graph <- merge(umap_graph, cluster_assignment, by.x="from", by.y="Row.names")
umap_graph <- umap_graph[umap_graph$cluster %in% c(4,13,17,19),]

write_csv(umap_graph, "/home/degan/ip_proteomics/inputs/umap_graph_ip.csv")


# plot expression of specific genes
umap_mat <- cds@reduce_dim_aux@listData[["UMAP"]]@listData[["model"]]@listData[["umap_model"]][["embedding"]]
umap_mat <- merge(umap_mat, cluster_assignment, by.x=0, by.y="Row.names")
umap_mat <- umap_mat[umap_mat$cluster == 13,]

rownames(umap_mat) <- NULL
umap_mat <- column_to_rownames(umap_mat, "Gene.names")
umap_dist <- as.matrix(dist(umap_mat, method='euclidean'))
umap_dist <- 1 - umap_dist
umap_dist <- umap_dist[!rownames(umap_dist) == "",]

colfunc <- colorRampPalette(c("white","#ffe6e6","#ffd4d4","#FFD5D5","#FFAAAA","#FF8080","#FF0000", "#ff0808"))

breaksList = seq(-1, 1, by = 0.001)
pheatmap(umap_dist, show_colnames = F, fontsize_row = 6.5,  color = colfunc(length(breaksList)),
         breaks = breaksList, filename = "/home/degan/ip_proteomics/figures/UMAP/dist_to_pd1.pdf", height = 5.5, width = 5.5,treeheight_row = 30,
         treeheight_col = 10)


## pathway plot ####
pd1_gprofiler <- read_xlsx("/home/degan/ip_proteomics/inputs/pd1_gprofiler.xlsx",skip = 1, trim_ws = TRUE)
pd1_gprofiler <- column_to_rownames(pd1_gprofiler,"Pathway")
pd1_gprofiler <- pd1_gprofiler %>% mutate_all(as.numeric)
pd1_gprofiler <- -log10(pd1_gprofiler)

breaksList = seq(0, 10, by = 0.01)
pheatmap::pheatmap(as.matrix(pd1_gprofiler),cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(brewer.pal(n = 5, name =
                                                         "Blues"))(length(breaksList)), breaks = breaksList,
                   filename = "/home/degan/ip_proteomics/figures/UMAP/GO_term_each_cluster.pdf", height = 4.5, width = 5.5)
                   