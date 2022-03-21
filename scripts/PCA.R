## library ####
library(prcomp)
library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)
library(PCA)
library(ggpubr)
library(reshape2)
library(gridExtra)

#### Note : stimulation: aCD3 aCD28 cell:bead 1:1 ## 

set.seed(16351893)
## Read in imputed data matrix and anntotations ####
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

## add imputed pd1 values to annotation ####
annotation_df <- cbind(annotation_df, t(imputed_data["Q15116", , drop = FALSE]))
names(annotation_df)[6] <- "pd1_imputed"

## Running PCA ####
pca <- prcomp(t(imputed_data), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)
#PCAloadings <- PCAloadings %>% arrange(desc(abs(PC6)))

## adding in ctrl vs stim variable ####
PCAvalues$antibody <- annotation_df$antibody 
PCAvalues$condition <- annotation_df$condition 
PCAvalues$time <- annotation_df$time
PCAvalues$stim <- annotation_df$stimulation
PCAvalues$replicate <- annotation_df$replicate
PCAvalues$pdcd1 <- log2(annotation_df$pd1_imputed)

## pca for each antibody condition ####
pdf("/home/degan/ip_proteomics/figures/pca_stim&time.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = antibody, group = condition)) +
  geom_point(aes(shape = condition), size = 2) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()


## running PCA for each antibody separately ####

### Nivo ####
nivo_annotation <- annotation_df[which(annotation_df$antibody %in% c("Niv")),]
nivo_imputed <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
nivo_imputed <- nivo_imputed[,which(colnames(nivo_imputed) %in% rownames(nivo_annotation))]

## before crap removal ####
pca <- prcomp(t(nivo_imputed), scale. = T, center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- nivo_annotation$condition
PCAvalues$replicate <- nivo_annotation$replicate
PCAvalues$stim <- nivo_annotation$stimulation
PCAvalues$time <- nivo_annotation$time

pdf("/home/degan/ip_proteomics/figures/pca_con&time_nivo.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = time)) +
  geom_point(aes(shape = time), size = 2) + scale_color_brewer(type = "qual")  +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

### Pem ####

pem_annotation <- annotation_df[which(annotation_df$antibody %in% c("Pem")),]
pem_imputed <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
pem_imputed <- pem_imputed[,which(colnames(pem_imputed) %in% rownames(pem_annotation))]

## before crap removal ####
pca <- prcomp(t(pem_imputed), scale. = T, center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- pem_annotation$condition
PCAvalues$replicate <- pem_annotation$replicate
PCAvalues$stim <- pem_annotation$stimulation
PCAvalues$time <- pem_annotation$time

pdf("/home/degan/ip_proteomics/figures/pca_con&time_pem.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = time)) +
  geom_point(aes(shape = time), size = 2) + scale_color_brewer(type = "qual", palette = "Set2") +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()



### Csab ####

csab_annotation <- annotation_df[which(annotation_df$antibody %in% c("Csa")),]
csab_imputed <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
csab_imputed <- csab_imputed[,which(colnames(csab_imputed) %in% rownames(csab_annotation))]

## before crap removal ####
pca <- prcomp(t(csab_imputed), scale. = T, center = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- csab_annotation$condition
PCAvalues$replicate <- csab_annotation$replicate
PCAvalues$stim <- csab_annotation$stimulation
PCAvalues$time <- csab_annotation$time

pdf("/home/degan/ip_proteomics/figures/pca_con&time_csab.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = time)) +
  geom_point(aes(shape = time), size = 2) + scale_color_brewer(type = "qual", palette = "Set1") +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()






