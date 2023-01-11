library(ggplot2)
library(reshape2)

## Read in imputed data matrix and anntotations ####
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")
pgroups.meta <- readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds")

## 1 calculate medians for each protein ####
medians_proteins <- apply(imputed_data, 1, median)
sample_deviations <- sweep(imputed_data, 2, medians_proteins)

## split in two for visualization purposes #####
cohort_1 <- sample_deviations[,c(1:120)]
cohort_2 <- sample_deviations[,c(121:240)]

cohort_1 <- melt(cohort_1)
cohort_2 <- melt(cohort_2)

pdf("/home/degan/ip_proteomics/figures/1_rle_plots.pdf", width = 6, height = 4)
ggplot(cohort_1, aes(x=variable, y=value)) + geom_boxplot(outlier.shape = NA) +
  theme_bw() + theme(axis.text.x = element_blank()) + xlab("Sample") + ylab("Relative Abundance") + 
  geom_hline(yintercept = 0, linetype = "dashed", color="Red")
dev.off()

pdf("/home/degan/ip_proteomics/figures/2_rle_plots.pdf", width = 6, height = 4)
ggplot(cohort_2, aes(x=variable, y=value)) + geom_boxplot(outlier.shape = NA) +
  theme_bw() + theme(axis.text.x = element_blank()) + xlab("Sample") + ylab("Relative Abundance") +
  geom_hline(yintercept = 0, linetype = "dashed", color="Red")
dev.off()
