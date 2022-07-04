library("leiden")
library(WGCNA)
library(igraph)
library("RColorBrewer")
library(factoextra)

## read in protein data ####
################################################################################

protein_data <- readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")
protein_data_filtered <- protein_data[which(rownames(protein_data) %in% diff_exp_filter$rn),]

## perform PCA on proteins
################################################################################

pca <- prcomp(protein_data, center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

ggplot(PCAvalues, aes(x = PC1, y = PC2)) +
  geom_point(size = 1) + theme_bw() + ggplot2::labs(x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),
                                                    y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))

fviz_eig(pca)


## WGCNA to create adj matrix ####
################################################################################

## choosing a soft-thresholding power value ####
powers = c(c(1:10), seq(from = 12, to=25, by=2))
sft = pickSoftThreshold(t(PCAvalues), powerVector = powers, verbose = 5)
#Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## picked soft power of 3, which is the lowest power for which the scale-free 
## topology fit index reaches 0.90
softPower = 7
adjacency = adjacency(t(PCAvalues), power = softPower)
adjacency <- as.matrix(adjacency)

## Leiden clustering algorithm ####
################################################################################

graph_object <- graph_from_adjacency_matrix(adjacency, mode = "directed", weighted = T,
                                            diag = FALSE)

partition <- leiden(adjacency, resolution_parameter = 0.8, n_iterations = -1)
table(partition)

node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols, vertex.label=NA, size = 5, arrow.size=0.2)

member_assignment <- data_frame(cluster = partition, row.names = rownames(PCAvalues))


