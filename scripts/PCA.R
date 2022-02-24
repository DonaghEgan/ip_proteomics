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

set.seed(16351893)
## Read in imputed data matrix and anntotations ####
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

## add imputed pd1 values to annotation ####
annotation_df <- cbind(annotation_df, t(imputed_data["Q15116", , drop = FALSE]))
names(annotation_df)[6] <- "pd1_imputed"

## Running PCA ####
pca <- prcomp(t(imputed_data))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)
#PCAloadings <- PCAloadings %>% arrange(desc(abs(PC6)))

## adding in ctrl vs stim variable ####
PCAvalues$antibody <- annotation_df$antibody 
PCAvalues$condition <- annotation_df$condition 
PCAvalues$time <- annotation_df$time
PCAvalues$batch <- annotation_df$batch
PCAvalues$pdcd1 <- annotation_df$pd1_imputed

## pca for each antibody condition ####
pdf("/home/degan/ip_proteomics/figures/pca_condtion&antibody.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = antibody, group = condition)) +
  geom_point(aes(shape = condition), size = 2) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

## removing stimulation experiments are rerunning ####
annotation_df$stimulation <- ifelse(grepl("Stim",rownames(annotation_df)), "Stimulated", "NotStim")
annotation_nostim <- annotation_df[which(annotation_df$stimulation == "NotStim"),]

imputed_nostim <- imputed_data[,which(colnames(imputed_data) %in% rownames(annotation_nostim))]

## re-running PCA ####
pca <- prcomp(t(imputed_nostim))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)
#PCAloadings <- PCAloadings %>% arrange(desc(abs(PC6)))

## adding in ctrl vs stim variable ####
PCAvalues$antibody <- annotation_nostim$antibody 
PCAvalues$condition <- annotation_nostim$condition 
PCAvalues$time <- annotation_nostim$time
PCAvalues$replicate <- annotation_nostim$replicate
PCAvalues$pdcd1 <- annotation_nostim$pd1_imputed

pdf("/home/degan/ip_proteomics/figures/pca_condtion&antibody_nostim_crapremove.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC6, colour = antibody, group = condition)) +
  geom_point(aes(shape = condition), size = 2) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

ggplot(PCAvalues, aes(x = PC1, y = PC6, colour = condition, group = replicate)) +
  geom_point(aes(shape=replicate), size = 2) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

pdf("/home/degan/ip_proteomics/figures/pca_pdcd1exp_nostim.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = pdcd1)) +
  geom_point(size = 2) + scale_color_viridis(discrete = F) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

pdf("/home/degan/ip_proteomics/figures/pca_batch_nostim.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = batch)) +
  geom_point(size = 2) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()


## running PCA for each antibody separately ####

### Nivo ####
nivo_annotation <- annotation_nostim[which(annotation_nostim$antibody %in% c("ConNivo", "PD1Nivo")),]
nivo_imputed <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
nivo_imputed <- nivo_imputed[,which(colnames(nivo_imputed) %in% rownames(nivo_annotation))]
imputed_nostim <- imputed_nostim[which(rownames(imputed_nostim) %in% gene_name_matrix$Protein.IDs),]
nivo_imputed_crap <- imputed_nostim [,which(colnames(imputed_nostim) %in% rownames(nivo_annotation))]

## before crap removal ####
pca <- prcomp(t(nivo_imputed))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- nivo_annotation$condition
PCAvalues$replicate <- nivo_annotation$replicate

pca_before_crap <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = replicate)) +
  geom_point(aes(shape = replicate), size = 2) + scale_color_brewer(type = "qual")  +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

## after crap removal ####
pca <- prcomp(t(nivo_imputed_crap))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- nivo_annotation$condition
PCAvalues$replicate <- nivo_annotation$replicate

pca_after_crap <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = replicate)) +
  geom_point(aes(shape = replicate), size = 2) + scale_color_brewer(type = "qual") +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

grid.arrange(pca_before_crap, pca_after_crap, ncol=2)


### Pem ####

pem_annotation <- annotation_nostim[which(annotation_nostim$antibody %in% c("ConPem", "PD1Pem")),]
pem_imputed <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
pem_imputed <- pem_imputed[,which(colnames(pem_imputed) %in% rownames(pem_annotation))]
imputed_nostim <- imputed_nostim[which(rownames(imputed_nostim) %in% gene_name_matrix$Protein.IDs),]
pem_imputed_crap <- imputed_nostim[,which(colnames(imputed_nostim) %in% rownames(pem_annotation))]

## before crap removal ####
pca <- prcomp(t(pem_imputed))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- pem_annotation$condition
PCAvalues$replicate <- pem_annotation$replicate


pca_before_crap <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = replicate)) +
  geom_point(aes(shape = replicate), size = 2) + scale_color_brewer(type = "qual") +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

## after crap removal ####
pca <- prcomp(t(pem_imputed_crap))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- pem_annotation$condition
PCAvalues$replicate <- pem_annotation$replicate

pca_after_crap <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = replicate)) +
  geom_point(aes(shape = replicate),size = 2) + scale_color_brewer(type = "qual") +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

grid.arrange(pca_before_crap, pca_after_crap, ncol=2)

### Csab ####

csab_annotation <- annotation_nostim[which(annotation_nostim$antibody %in% c("ConCsab", "PD1Csab")),]
csab_imputed <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
csab_imputed <- csab_imputed[,which(colnames(csab_imputed) %in% rownames(csab_annotation))]
imputed_nostim <- imputed_nostim[which(rownames(imputed_nostim) %in% gene_name_matrix$Protein.IDs),]
csab_imputed_crap <- imputed_nostim[,which(colnames(imputed_nostim) %in% rownames(csab_annotation))]

## before crap removal ####
pca <- prcomp(t(csab_imputed))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- csab_annotation$condition
PCAvalues$replicate <- csab_annotation$replicate

pca_before_crap <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = replicate)) +
  geom_point(aes(shape = replicate), size = 2) + scale_color_brewer(type = "qual") +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

## after crap removal ####
pca <- prcomp(t(csab_imputed_crap))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues$condition <- csab_annotation$condition
PCAvalues$replicate <- csab_annotation$replicate

pca_after_crap <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, group = replicate)) +
  geom_point(aes(shape = replicate), size = 2) + scale_color_brewer(type = "qual") +
  theme_bw() +
  ggplot2::labs(x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),
                y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

grid.arrange(pca_before_crap, pca_after_crap, ncol=2)




