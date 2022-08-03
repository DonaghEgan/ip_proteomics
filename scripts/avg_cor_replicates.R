################################################################################
# Donagh Egan
# POIAZ 
# Date: 22nd June 2022
# 

# Description: This script calculates the correlation of replicates (4 biological x 2 technical)
# within each experimental condition
################################################################################

## Library #### 
################################################################################

## Load in protein data ####
################################################################################
protein_data <- readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")

## read in exp design - created by Vadim ####
experimental_design <- read.csv(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/experimental_design1.csv", sep = "\t", 
                                colClasses=c("character","character"), header = F)

experimental_design$batch = sapply(strsplit(colnames(protein_data), "_"), "[", 4)
colnames(experimental_design) <- c("label","condition", "batch")

##remove stimulation from exp design
experimental_design$stimulation <- ifelse(grepl("Stim",experimental_design$label), "Stimulated", "NotStim")
experimental_design <- experimental_design[which(experimental_design$stimulation == "NotStim"),]

## calculate pairwise correlation between all samples ####
################################################################################

cor_protein_data <- cor(protein_data)
cor_protein_data <- merge(cor_protein_data, experimental_design, by.x=0, by.y="label")
cor_protein_data <- column_to_rownames(cor_protein_data, "Row.names")

## subset correlations according to conditions ####
################################################################################
list_conditions <- list()
for (i in unique(cor_protein_data$condition)) {
  df_temp <- cor_protein_data[cor_protein_data$condition == i, , drop=F]
  df_temp <- df_temp[, which(colnames(df_temp) %in% rownames(df_temp))]
  diag(df_temp) <- NA
  list_conditions[[i]] <- df_temp
  
}

## find mean corrlation/conditon and add to list 
################################################################################
cor_values <- list()
for (i in names(list_conditions)) {
  values <- list_conditions[[i]][,1]
  values <- mean(values, na.rm=T)
  sd <- sd(values)
  cor_values[[i]] <- values
}

## Format and plot results ####
################################################################################

cor_df <- data.frame(avg_cor = matrix(unlist(cor_values), nrow=length(cor_values), byrow=TRUE),
                     row.names = names(cor_values))
cor_df$condition <- rownames(cor_df)
nb.cols <- 30
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

pdf("/home/degan/ip_proteomics/figures/pre_process/avg_corr_condition.pdf", width = 5, height = 4)
ggplot(data=cor_df, aes(x=condition, y=avg_cor, fill = condition)) +
  geom_bar(stat="identity", width = 0.7, fill=pal_npg("nrc")(10)[1]) +  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(aes(yintercept = mean(avg_cor)), linetype = "dashed") + ylim(0,1) +
  ylab("Average PCC")+ xlab("Condition") + theme(legend.position = "none")
dev.off()

## correlation according to batch ####
################################################################################

## subset correlations according to batch ####
################################################################################
list_batch <- list()
for (i in unique(cor_protein_data$batch)) {
  df_temp <- cor_protein_data[cor_protein_data$batch == i, , drop=F]
  df_temp <- df_temp[, which(colnames(df_temp) %in% rownames(df_temp))]
  diag(df_temp) <- NA
  list_batch[[i]] <- df_temp
  
}

cor_batch <- list()
for (i in names(list_batch)) {
  values <- list_batch[[i]][,1]
  values <- mean(values, na.rm=T)
  sd <- sd(values)
  cor_batch[[i]] <- values
}

cor_batch_df <- data.frame(avg_cor = matrix(unlist(cor_batch), nrow=length(cor_batch), byrow=TRUE),
                     row.names = names(cor_batch))

res <- wilcox.test(cor_df$avg_cor, cor_batch_df$avg_cor)
pdf("/home/degan/ip_proteomics/figures/pre_process/cor_batch_replicates.pdf", width = 3, height = 3.5)
boxplot(cor_df$avg_cor, cor_batch_df$avg_cor, col= pal_npg("nrc",alpha = 0.7)(10)[3:5],
        ylim = c(0.5, 1), ylab = "Pearson's r", cex.axis=0.7, cex.lab=0.8)
text(x=1.5, y=0.95,
     paste("Mann Whitney: ", round(res$p.value,digits = 4), sep=""), pos=3, cex=0.7)
dev.off()

## random correlations ####
################################################################################
random <-  sample(cor_protein_data, replace = T, 1000)

res_batch_random <- wilcox.test(colMeans(random[sapply(random, is.numeric)]), cor_batch_df$avg_cor)
res_random_rep <-  wilcox.test(colMeans(random[sapply(random, is.numeric)]), cor_df$avg_cor)
res_batch_rep <- wilcox.test(cor_df$avg_cor, cor_batch_df$avg_cor)

pdf("/home/degan/ip_proteomics/figures/pre_process/cor_batch_rep_random.pdf", width = 3, height = 3.5)
boxplot(cor_df$avg_cor, cor_batch_df$avg_cor,colMeans(random[sapply(random, is.numeric)]), col= pal_npg("nrc",alpha = 0.7)(10)[1:5],
        ylim = c(0.5, 1), ylab = "Mean Pearson's r", cex.axis=0.7, cex.lab=0.8)
text(x=1, y=0.93, "***", pos=3, cex=1)
dev.off()





