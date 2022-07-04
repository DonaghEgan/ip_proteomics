################################################################################
# Donagh Egan
# POIAZ 
# Date: 14 March 2022
# 

# Description: Preforms WGCNA on ip proteomics data using proteins identified from Limma
#              as differentially bound across PD1 time-points
################################################################################

## Library #### 
################################################################################

library(stringi)
library(WGCNA)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(mgsa)
library(clusterProfiler)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

## Load in filtered differentially bound data and processed counts #####
diff_bound_filtered <- readRDS("/home/degan/ip_proteomics/inputs/diff_bound_filter.Rds")
processed_ip <- readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")
experimental_design <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

## gene list for each comparison made ####
gene_list <- list(PD1vsCtr = diff_bound_filtered$rn[diff_bound_filtered$type=="PD1vsCtr"],
                  "PD1-time" = diff_bound_filtered$rn[diff_bound_filtered$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")],
                  "PD1vsCtr-time" = diff_bound_filtered$rn[diff_bound_filtered$type %in% c("PD1_0vsCon_0", "PD1_20vsCon_20", 
                                                                                       "PD1_24vsCon_24",   "PD1_5vsCon_5")])

## top differentially bound proteins across time within PD1 ####
top_bound_pd1time <- diff_bound_filtered[(diff_bound_filtered$rn %in% gene_list[["PD1-time"]] &
                                        diff_bound_filtered$type %in% c("ST5vNsT0", "ST20vST5", "ST24vST20")),] %>% 
  arrange(desc(t))


processed_bound_pd1time <- processed_ip[which(rownames(processed_ip) %in% top_bound_pd1time$rn),]

## choosing a soft-thresholding power value ####
powers = c(c(1:10), seq(from = 12, to=25, by=2))
sft = pickSoftThreshold(t(processed_bound_pd1time), powerVector = powers, verbose = 5)
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
softPower = 6
adjacency = adjacency(t(processed_bound_pd1time), power = softPower)

## Turn adjacency into topological overlap matrix ####
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

## Call the hierarchical clustering function ####
geneTree = hclust(as.dist(dissTOM), method = "complete")

## Plot the resulting clustering tree (dendrogram) ####
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

## Module identification using dynamic tree cut ####
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 5)
table(dynamicMods)

## Convert numeric lables into colors ####
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

## Plot the dendrogram and colors underneath ####
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

## Calculate eigengenes ####
MEList = moduleEigengenes(t(diff_exp_counts_pd1), colors = dynamicColors)
MEs = MEList$eigengenes
## Calculate dissimilarity of module eigengenes ####
MEDiss = 1-cor(MEs)
## Cluster module eigengenes ####
METree = hclust(as.dist(MEDiss), method = "complete")
## Plot the result ####
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.5
## Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

## call an automatic merging function
merge = mergeCloseModules(t(diff_exp_counts_pd1), dynamicColors, cutHeight = MEDissThres, verbose = 3)

## The merged module colors
mergedColors = merge$colors

## Eigengenes of the new merged modules
mergedMEs = merge$newMEs

## plotting gene dendogram with original and merged module colors underneath ####
sizeGrWindow(12, 9)
pdf(file = "/home/degan/ip_proteomics/figures/WGCNA/WGCNA_geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## Rename to moduleColors ####
moduleColors = mergedColors

## plotting module heatmap - all proteins #### 
sizeGrWindow(9,9)
TOMplot(dissTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

tom_annotation <- data.frame("module"=as.factor(moduleColors))
tom_annotation$ind <- rownames(tom_annotation)
tom_annotation <- data.frame(tom_annotation[geneTree[["order"]],])
rownames(tom_annotation) <- paste("X", rownames(tom_annotation), sep="") 
tom_annotation$ind <- NULL

mycolors <- c(unique(as.character(tom_annotation$module)))
names(mycolors) <-  unique(tom_annotation$module)
mycolors <- list(module = mycolors)

TOM_hm <- data.frame(dissTOM)
pdf("/home/degan/ip_proteomics/figures/WGCNA/module_heatmap.pdf", height = 4, width = 4)
pheatmap::pheatmap(TOM_hm, clustering_distance_cols = "euclidean", annotation = tom_annotation,
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                   clustering_method = "complete", annotation_colors = mycolors,
                   show_rownames = F, show_colnames = F, treeheight_row = 12,
                   treeheight_col = 15)
dev.off()



## Construct numerical labels corresponding to the colors ####
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

## Recalculate MEs with color labels ####
nSamples = nrow(experimental_design)

MEs0 = moduleEigengenes(t(processed_bound_pd1time), moduleColors)$eigengenes
MEs = orderMEs(MEs0)
experimental_design$time <- as.numeric(stri_replace_all_regex(experimental_design$time, 
                                                              pattern = c('min', 'hr'), 
                                                              replacement = "", vectorize = F))

moduleTraitCor = cor(MEs, experimental_design[,c(3),drop=F], use = "p", method = "pearson")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(experimental_design[,c(3),drop=F]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = rev(greenWhiteRed(50)),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## protein relationship to trait and important modules ####

## Define variable time 
time = as.data.frame(experimental_design$time)
names(time) = "time"
## names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(t(processed_bound_pd1time), MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(t(processed_bound_pd1time), time, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(time), sep="")
names(GSPvalue) = paste("p.GS.", names(time), sep="")

## plotting signficance of a genes correlation with time vs significance of a genes
## correlation with black module 
module = "tan"
column = match(module, modNames)
moduleGenes = moduleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

## Get the corresponding Locis Link IDs ##
hs <- org.Hs.eg.db
uniprot_id <- sapply(strsplit(rownames(processed_bound_pd1time), ";"), "[", 1)

## find module-specific proteins ####
intModules <- c(gsub("ME", "", names(MEs))) 
modGeneslist = list()
for (module in intModules) {
  modGenes = (moduleColors==module)
  modLLIDs = uniprot_id[modGenes]
  modGeneslist[[module]] = modLLIDs
}

## formatting GO annotations ####
## read in file
GO <- readGAF("/home/degan/phospho_proteomics/inputs/goa_human.gaf")

## extract relevant info
mapping.index <-  GO@itemName2ItemIndex
ID.annotations <- itemAnnotations(GO)
GO.sets <- GO@sets
GO.annotation <- setAnnotations(GO)

## create df - GOID and ID index
GO.df <- data.frame("GOID" = rep(names(GO.sets), sapply(GO.sets, length)),
                    "ID.index" = unlist(GO.sets),  row.names = NULL)

## Remove category 'all', and add column with Uniprot ids
GO.annotation <- GO.annotation[GO.annotation[,"term"] != "all", ]
GO.annotation[,"GOID"] <- rownames(GO.annotation)
# GO.df
GO.df <- GO.df[GO.df[,"GOID"] != "all", ]
GO.df[,"UNIPROTKB"] <- names(mapping.index [GO.df[,"ID.index"] ])

## background proteins ####
proteins.all <- sapply(strsplit(rownames(processed_ip), ";"), "[", 1)

res.GO.ora_neg <- enricher(
  gene=c(modGeneslist[["grey"]]),
  pvalueCutoff = 1,
  pAdjustMethod = "bonferroni",
  universe=proteins.all,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 1,
  TERM2GENE = GO.df[ ,c("GOID","UNIPROTKB")],
  TERM2NAME = GO.annotation[ ,c("GOID", "term")]
)

res.GO.ora_pos <- enricher(
  gene=c(modGeneslist[["pink"]], modGeneslist[["grey"]]),
  pvalueCutoff = 1,
  pAdjustMethod = "bonferroni",
  universe=proteins.all,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 1,
  TERM2GENE = GO.df[ ,c("GOID","UNIPROTKB")],
  TERM2NAME = GO.annotation[ ,c("GOID", "term")]
)

## generate plots dotplot ####
library(enrichplot)

pdf("/home/degan/ip_proteomics/figures/WGCNA/mod_neg_cor.pdf", height = 4, width = 6)
dotplot(res.GO.ora_neg, showCategory=9,font.size = 8)
dev.off()

pdf("/home/degan/ip_proteomics/figures/WGCNA/mod_pos_cor.pdf", height = 4, width = 6)
dotplot(res.GO.ora_pos, showCategory=9,font.size = 8)
dev.off()

cols_cor <- c("MEpink","MEgrey", "MEtan", "MEred")
experimental_design <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")

pdf("/home/degan/ip_proteomics/figures/WGCNA/module_heatmap_time.pdf", height = 5, width = 5)
pheatmap::pheatmap(t(MEs[,cols_cor]), scale = "row", annotation = experimental_design[,c(2,3), drop=FALSE],
         show_colnames = F, clustering_distance_cols = "correlation", treeheight_row = 4,
         treeheight_col = 5)
dev.off()








