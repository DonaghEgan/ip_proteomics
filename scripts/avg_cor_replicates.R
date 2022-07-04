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


## calculate pairwise correlation between all samples ####
################################################################################

cor_protein_data <- cor(protein_data)
cor_protein_data <- merge(cor_protein_data, experimental_design, by.x=0, by.y = "label")
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
  geom_bar(stat="identity", width = 0.7) + scale_fill_manual(values = mycolors) +  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(aes(yintercept = mean(avg_cor)), linetype = "dashed") + ylim(0,1) +
  ylab("Average PCC")+ xlab("Condition") + theme(legend.position = "none")
dev.off()

