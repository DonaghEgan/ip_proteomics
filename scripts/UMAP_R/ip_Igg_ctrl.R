library(readxl)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
library(tidyverse)
library(tidyr)

DIA_PD_1_sig <- read_excel("ip_proteomics/inputs/2022 07 08 exp83 DIA PD-1 sig.xlsx", 
                           sheet = "comp")

cluster_assignment <- cluster_assignment %>% mutate(Gene.names = strsplit(as.character(Gene.names), ";")) %>% tidyr::unnest(Gene.names)

colnames(DIA_PD_1_sig)[1] <- "Common Proteins"
DIA_PD_1_sig <- na.omit(DIA_PD_1_sig[ ,"Common Proteins"])

GENEIDS_H <- AnnotationDbi::select(hs,keys = DIA_PD_1_sig$`Common Proteins`, columns = c("UNIPROT", "SYMBOL"), keytype = "UNIPROT")
DIA_PD_1_sig <- merge(DIA_PD_1_sig, GENEIDS_H, by.x = "Common Proteins", by.y = "UNIPROT")

intersect(cluster_assignment$Gene.names, DIA_PD_1_sig$SYMBOL)

cluster_assignment_sub <- cluster_assignment[cluster_assignment$Gene.names %in% DIA_PD_1_sig$SYMBOL,]

cluster_overlap_freq <- data.frame(table(cluster_assignment_sub$cluster),
                                   cluster_size = table(cluster_assignment$cluster))

cluster_overlap_freq$percent_overlap <- (cluster_overlap_freq$Freq / cluster_overlap_freq$cluster_size.Freq) * 100
