top_proteins <- readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")

gene_name_matrix <- pgroups[which(pgroups$Protein.IDs %in% rownames(top_proteins)),]
top_proteins <- cbind(top_proteins, gene_name_matrix[, "Gene.names", drop=FALSE])
rownames(top_proteins) <- NULL
top_proteins <- column_to_rownames(top_proteins, "Gene.names")

## top proteins in PD1 samples only ####
pd1_only <- annotation_nostim[which(substr(annotation_nostim$antibody, 1, 3) == "PD1"),]
top_proteins_pd1 <- top_proteins[,which(colnames(top_proteins) %in% rownames(pd1_only))]

cor_proteins <- cor(t(top_proteins_pd1), method = "spearman")

pdf("/home/degan/ip_proteomics/figures/PDCD1_module_heatmap.pdf", height = 5, width = 5)
protein_heatmap <- ComplexHeatmap::Heatmap(cor_proteins,name = "Corr", show_row_names = T, 
                        show_column_names = FALSE, row_names_gp = gpar(fontsize = 4))
dev.off()

protein_modules <- draw(protein_heatmap)

PDCD1_module <- c("LONP1", "CAAP1", "TWISTNB", "SUV39H1", "SNIP1", "GLTSCR2",
                  "SFSWAP", "EXOSC5", "SRSFB", "SREK1IP1", "PDCD1", "AHNAK",
                  "SATB1", "MSANTD4", "DDX42", "SP1", "WHSC1L1", "TERF1", "SS18",
                  "NSMCE2", "WDR55", "PML", "SCAF1", "PUF60", "CDC40", "RRP1B",
                  "MAX", "USF1", "DR1", "RREB1", "NACA", "DIAPH1", "PTPN11", "DNAJC2",
                  "ARRB2", "RRN3:RRN3P3", "INTS3", "OR10AG1", "SUN1", "BCL11A",
                  "URB2", "SMARCA2", "HIST1H2BL;HIST1H2BA", "MDC1", "CENPV",
                  "CASP8AP2", "C19orf47", "SSSCA1")

top_proteins_pd1trans <- data.frame(t(top_proteins_pd1))

PDCD1_module_time <- top_proteins_pd1trans[,which(colnames(top_proteins_pd1trans) %in% PDCD1_module)]
PDCD1_module_time$Time <- sapply(strsplit(rownames(PDCD1_module_time), "_"), "[", 2)

PDCD1_module_time_mean <- aggregate(. ~ Time, PDCD1_module_time, median)
PDCD1_module_time_mean <- melt(PDCD1_module_time_mean, id.vars = "Time")
PDCD1_module_time_mean$Time <- gsub("min", "", PDCD1_module_time_mean$Time) 
PDCD1_module_time_mean$Time <- gsub("hr", "", PDCD1_module_time_mean$Time) 
PDCD1_module_time_mean$Time <- as.numeric(PDCD1_module_time_mean$Time)
PDCD1_module_time_mean$Time <- ifelse(PDCD1_module_time_mean$Time == 24, 
                                      PDCD1_module_time_mean$Time*60, 
                                      PDCD1_module_time_mean$Time*1)

PDCD1_module_time_mean <- column_to_rownames(PDCD1_module_time_mean, "Time")
PDCD1_module_time_mean <- data.frame(t(PDCD1_module_time_mean))

pheatmap(t(PDCD1_module_time[,1:45]), annotation_col = PDCD1_module_time[,46, drop=F],
         cluster_cols = T, cluster_rows = T)







