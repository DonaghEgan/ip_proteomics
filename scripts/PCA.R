## library ####
library(prcomp)
library(ggplot2)
library(dplyr)
library(viridis)

protein_matrix <- read.delim("~/ip_proteomics/inputs/matrix13_Perseus.txt", header=TRUE,)
protein_matrix = protein_matrix[!duplicated(protein_matrix$T..Gene.names),] ## remove duplicate genes 
rownames(protein_matrix) <- protein_matrix$T..Gene.names

## LFQ intensities only ####
protein_matrix <- protein_matrix[,1:180]

## Running PCA ####
pca <- prcomp(t(protein_matrix))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)
PCAloadings <- PCAloadings %>% arrange(desc(abs(PC6)))

## adding in ctrl vs stim variable ####
PCAvalues$condition <- ifelse(grepl("Con",rownames(PCAvalues)), "Control", "Antibody")

# Plot
pdf("/home/degan/ip_proteomics/figures/PCA&Loadings.pdf", width = 4, height = 3.5)
ggplot(PCAvalues, aes(x = PC1, y = PC6, colour = condition)) +
  geom_segment(data = PCAloadings[1:4,], aes(x = 0, y = 0, xend = (PC1*80),
                                       yend = (PC6*80)), arrow = arrow(length = unit(.25, "picas")),
               color = "black") +
  geom_point(size = 1) + scale_color_viridis(discrete = TRUE) +
  ggrepel::geom_label_repel(data=PCAloadings[1:4,],aes(x= (PC1*80), y = (PC6*80),
                                                 label = Variables),
                            size = 1.7,
                            fontface = "bold",
                            show.legend = FALSE,
                            inherit.aes = FALSE,
  ) + theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y = paste0("PC6: ", round(percentVar[2] * 100), "% variance"))
dev.off()


