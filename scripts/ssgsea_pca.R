library(GSVA)
library(tidyverse)
library(viridis)


## Read in imputed data matrix and anntotations ####
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

pgroups <- read.table(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/proteinGroups_PPI_Martina_Jan.txt", 
                      header = T, sep = "\t",quote='',stringsAsFactors = FALSE,comment.char="")

pgroups.meta <- pgroups[,c("Protein.IDs","Gene.names")]


## Setting up pathways ####
msigdbr_df <- msigdbr(species = "human", category = "C2", subcategory = "CP")
pathways = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

## adding gene names imputed_data ####
imputed_data <- merge(imputed_data, pgroups.meta, by.x = 0, by.y = "Protein.IDs")
imputed_data <- imputed_data[!duplicated(imputed_data$Gene.names),]
imputed_data$Row.names <- NULL
rownames(imputed_data) <- NULL
imputed_data <- column_to_rownames(imputed_data, "Gene.names")

## ssgsea full data set ########################################################
################################################################################

## scroring 
set.seed(16351893)
genesets.gsva <- gsva(as.matrix(imputed_data), pathways, mx.diff=FALSE, verbose=TRUE, parallel.sz=1, method = "ssgsea") #ssgsea

## Running PCA ####
pca <- prcomp(t(genesets.gsva), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- annotation_df$condition 
PCAvalues$time <- annotation_df$time
PCAvalues$stim <- annotation_df$stimulation
PCAvalues$antibody <- annotation_df$antibody

pos_PC1 <- PCAloadings %>% arrange(desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1)  
neg_PC1 <- PCAloadings %>% arrange(-desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1) 
pos_PC2 <- PCAloadings %>% arrange(desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)
neg_PC2 <- PCAloadings %>% arrange(-desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)

df_list <- list(pos_PC1, neg_PC1, pos_PC2, neg_PC2)
PCAloadings_annot <- do.call(rbind, df_list)

## Plot ####
pdf("/home/degan/ip_proteomics/figures/pca_ssgsea_allanti.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = antibody)) +
  geom_segment(data = PCAloadings_annot, aes(x = 0, y = 0, xend = (PC1*10),
                                       yend = (PC2*10)), arrow = arrow(length = unit(.25, "picas")),
               color = "black") +
  geom_point(aes(shape=condition), size = 1.2) + scale_color_viridis(discrete = T) +
  ggrepel::geom_label_repel(data=PCAloadings_annot,aes(x= (PC1*10), y = (PC2*10),
                                                 label = Variables),
                            size = 2,
                            fontface = "bold",
                            show.legend = FALSE,
                            label.size = NA,
                            fill = NA,
                            inherit.aes = FALSE,
  ) +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw()
dev.off()

## ssgsea pem only #############################################################
################################################################################

pem_annotation <- annotation_df[which(annotation_df$antibody %in% c("Pem")),]
pem_imputed <- imputed_data[,which(colnames(imputed_data) %in% rownames(pem_annotation))]

## scroring ####
set.seed(16351893)
genesets.gsva_pem <- gsva(as.matrix(pem_imputed), pathways, mx.diff=FALSE, verbose=TRUE, parallel.sz=1, method = "ssgsea") #ssgsea

## Running PCA ####
pca <- prcomp(t(genesets.gsva_pem), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- pem_annotation$condition 
PCAvalues$time <- pem_annotation$time
PCAvalues$stim <- pem_annotation$stimulation
PCAvalues$antibody <- pem_annotation$antibody

pos_PC1 <- PCAloadings %>% arrange(desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1)  
neg_PC1 <- PCAloadings %>% arrange(-desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1) 
pos_PC2 <- PCAloadings %>% arrange(desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)
neg_PC2 <- PCAloadings %>% arrange(-desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)


df_list <- list(pos_PC1, neg_PC1, pos_PC2, neg_PC2)
PCAloadings_annot <- do.call(rbind, df_list)

pdf("/home/degan/ip_proteomics/figures/pca_ssgsea_pem_only.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition)) +
  geom_segment(data = PCAloadings_annot, aes(x = 0, y = 0, xend = (PC1*10),
                                             yend = (PC2*10)), arrow = arrow(length = unit(.25, "picas")),
               color = "black") +
  geom_point(size = 1.2) + scale_color_brewer(type = "qual", palette = "Set2") +
  ggrepel::geom_label_repel(data=PCAloadings_annot,aes(x= (PC1*10), y = (PC2*10),
                                                       label = Variables),
                            size = 2,
                            fontface = "bold",
                            show.legend = FALSE,
                            label.size = NA,
                            fill = NA,
                            inherit.aes = FALSE,
  ) +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw()
dev.off()

## ssgsea niv only #############################################################
################################################################################

niv_annotation <- annotation_df[which(annotation_df$antibody %in% c("Niv")),]
niv_imputed <- imputed_data[,which(colnames(imputed_data) %in% rownames(niv_annotation))]

## scroring ####
set.seed(16351893)
genesets.gsva_niv <- gsva(as.matrix(niv_imputed), pathways, mx.diff=FALSE, verbose=TRUE, parallel.sz=1, method = "ssgsea") #ssgsea

## Running PCA ####
pca <- prcomp(t(genesets.gsva_niv), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- niv_annotation$condition 
PCAvalues$time <- niv_annotation$time
PCAvalues$stim <- niv_annotation$stimulation
PCAvalues$antibody <- niv_annotation$antibody

pos_PC1 <- PCAloadings %>% arrange(desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1)  
neg_PC1 <- PCAloadings %>% arrange(-desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1) 
pos_PC2 <- PCAloadings %>% arrange(desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)
neg_PC2 <- PCAloadings %>% arrange(-desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)

df_list <- list(pos_PC1, neg_PC1, pos_PC2, neg_PC2)
PCAloadings_annot <- do.call(rbind, df_list)

pdf("/home/degan/ip_proteomics/figures/pca_ssgsea_niv_only.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition)) +
  geom_segment(data = PCAloadings_annot, aes(x = 0, y = 0, xend = (PC1*10),
                                             yend = (PC2*10)), arrow = arrow(length = unit(.25, "picas")),
               color = "black") +
  geom_point(aes(shape=time),size = 1.2) + scale_color_brewer(type = "qual") +
  ggrepel::geom_label_repel(data=PCAloadings_annot,aes(x= (PC1*10), y = (PC2*10),
                                                       label = Variables),
                            size = 2,
                            fontface = "bold",
                            show.legend = FALSE,
                            label.size = NA,
                            fill = NA,
                            inherit.aes = FALSE,
  ) +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw()
dev.off()


## ssgsea Csab only #############################################################
################################################################################

csab_annotation <- annotation_df[which(annotation_df$antibody %in% c("Csa")),]
csab_imputed <- imputed_data[,which(colnames(imputed_data) %in% rownames(csab_annotation))]

## scroring ####
set.seed(16351893)
genesets.gsva_csab <- gsva(as.matrix(csab_imputed), pathways, mx.diff=FALSE, verbose=TRUE, parallel.sz=1, method = "ssgsea") #ssgsea

## Running PCA ####
pca <- prcomp(t(genesets.gsva_csab), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- csab_annotation$condition 
PCAvalues$time <- csab_annotation$time
PCAvalues$stim <- csab_annotation$stimulation
PCAvalues$antibody <- csab_annotation$antibody

pos_PC1 <- PCAloadings %>% arrange(desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1)  
neg_PC1 <- PCAloadings %>% arrange(-desc(PC1)) %>% select(Variables,PC1, PC2) %>% slice(1) 
pos_PC2 <- PCAloadings %>% arrange(desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)
neg_PC2 <- PCAloadings %>% arrange(-desc(PC2)) %>% select(Variables,PC1, PC2) %>% slice(1)

df_list <- list(pos_PC1, neg_PC1, pos_PC2, neg_PC2)
PCAloadings_annot <- do.call(rbind, df_list)

pdf("/home/degan/ip_proteomics/figures/pca_ssgsea_csab_only.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition)) +
  geom_segment(data = PCAloadings_annot, aes(x = 0, y = 0, xend = (PC1*10),
                                             yend = (PC2*10)), arrow = arrow(length = unit(.25, "picas")),
               color = "black") +
  geom_point(aes(shape=time),size = 1.2) + scale_color_brewer(type = "qual", palette = "Set1") +
  ggrepel::geom_label_repel(data=PCAloadings_annot,aes(x= (PC1*10), y = (PC2*10),
                                                       label = Variables),
                            size = 2,
                            fontface = "bold",
                            show.legend = FALSE,
                            label.size = NA,
                            fill = NA,
                            inherit.aes = FALSE,
  ) +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw()
dev.off()














