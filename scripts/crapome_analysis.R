## Library ####
library(org.Hs.eg.db)
library(tidyverse)
library(AnnotationDbi)
library(ggsci)
library(scales)
library(stats)
library(viridis)
#

## read in file from crapome ####
jurkat_crap <- read.table(file = "/home/degan/ip_proteomics/inputs/1643709596_userCrapDB.xls", sep = "\t", header=TRUE)

## read in normalised and imputed data matrix ####
protein_data_imputed <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))

## removing stimulation experiments are rerunning ####
annotation_df$stimulation <- ifelse(grepl("Stim",rownames(annotation_df)), "Stimulated", "NotStim")
annotation_nostim <- annotation_df[which(annotation_df$stimulation == "NotStim"),]
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
imputed_nostim <- imputed_data[,which(colnames(imputed_data) %in% rownames(annotation_nostim))]

## extract gene names from original protein matrix ####
gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(imputed_nostim)),]
rownames(gene_name_matrix) <- NULL
gene_name_matrix <-  gene_name_matrix[ ! gene_name_matrix$Gene.names %in% jurkat_crap$GENE,]

## remove crapome values from imputed data matrix ####
protein_data_imputed <- protein_data_imputed[which(rownames(protein_data_imputed) %in% gene_name_matrix$Protein.IDs),]

## re-running PCA ####
pca <- prcomp(t(protein_data_imputed))
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

pdf("/home/degan/ip_proteomics/figures/pca_condtion&antibody.pdf", width = 5, height = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = antibody, group = condition)) +
  geom_point(aes(shape = condition), size = 2) + scale_color_viridis(discrete = T) +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

## clustering and annotation after crapome removal ####
col_antibody <- c(pal_npg("nrc")(6))
names(col_antibody) <- unique(annotation_nostim$antibody)

heat_annotation <- ComplexHeatmap::HeatmapAnnotation(antibody = annotation_nostim$antibody,
                                                     replicate = annotation_nostim$replicate,
                                                     time = annotation_nostim$time,
                                                     col = list(antibody = col_antibody),
                                                     show_annotation_name = FALSE)

## removing stimulation experiments and rerunning ####
annotation_df$stimulation <- ifelse(grepl("Stim",rownames(annotation_df)), "Stimulated", "NotStim")
annotation_nostim <- annotation_df[which(annotation_df$stimulation == "NotStim"),]
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
imputed_nostim <- imputed_data[,which(colnames(imputed_data) %in% rownames(annotation_nostim))]

mat_scaled = scale(imputed_nostim)
ComplexHeatmap::Heatmap(mat_scaled[1:100,], clustering_distance_columns = "euclidean",
                        show_row_names = FALSE,  top_annotation = heat_annotation, 
                        show_column_names = FALSE)



