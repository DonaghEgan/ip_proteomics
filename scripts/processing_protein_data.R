## Library ####
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(pheatmap)
library(tidyverse)
library(imputeLCMD)
library(vsn)
library(readr)

## Reading in file ####
pgroups <- read.table(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/proteinGroups_PPI_Martina_Jan.txt", 
                      header = T, sep = "\t",quote='',stringsAsFactors = FALSE,comment.char="")

## Removing contaminants and unnecessary characters ####
pgroups <- pgroups[pgroups$Reverse != "+" & pgroups$Potential.contaminant != "+", ]
pgroups <- pgroups[,grepl("*rotein.IDs*|*names*|LFQ*", colnames(pgroups))]

## meta data for transfering gene labels ####
pgroups.meta <- pgroups[,c("Protein.IDs","Gene.names")]
saveRDS(pgroups.meta, "/home/degan/ip_proteomics/inputs/pggroups.meta.Rds")

## read in exp design - created by Vadim ####
experimental_design <- read.csv(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/experimental_design1.csv", sep = "\t", 
                                colClasses=c("character","character"), header = F)
colnames(experimental_design) <- c("label","condition")

## setting rownames ####
rownames(pgroups) <- NULL
protein_data <- column_to_rownames(pgroups, "Protein.IDs")

## removing "LFQ" string and setting 0 values to NA ####ins
protein_data <- protein_data[,grep("LFQ", colnames(protein_data))]
protein_data[protein_data=="0"]<-NA
colnames(protein_data) <- gsub(x = colnames(protein_data), pattern = "LFQ\\.intensity\\.", replacement = "") 
colnames(protein_data) <- gsub('.{7}$', '', colnames(protein_data))

## creating annotation dataframe according to substrings ####
annotation_df <- data.frame(condition = substr(sapply(strsplit(colnames(protein_data), "_"), "[", 1), 1, 3), 
                            row.names =  colnames(protein_data),
                            antibody = substr(sapply(strsplit(colnames(protein_data), "_"), "[", 1), 4,6),
                            time = sapply(strsplit(colnames(protein_data), "_"), "[", 2),
                            batch = sapply(strsplit(colnames(protein_data), "_"), "[", 4),
                            stimulation = ifelse(grepl("Stim|_5min|_20",colnames(protein_data)),"stim","unstim"))

## Heatmap of missing values ####
missval <- protein_data[apply(protein_data, 1, function(x) any(is.na(x))), ] # select only proteins with NA values
missval <- ifelse(is.na(missval), 0, 1)

heat_annotation <- ComplexHeatmap::HeatmapAnnotation(time = annotation_df$time,
                                                     antibody =  annotation_df$antibody,
                                                     stim = annotation_df$stimulation,
                                                     batch = annotation_df$batch,
                                                     show_annotation_name = FALSE)

pdf("/home/degan/ip_proteomics/figures/missingval_pd1.pdf", height = 5, width = 10)
ComplexHeatmap::Heatmap(missval, col = c("white", "black"), column_names_side = "top", 
              show_row_names = FALSE,  top_annotation = heat_annotation, show_column_names = F, name = "Missing values pattern", 
              column_names_gp = gpar(fontsize = 3), show_column_dend = F, cluster_columns = T, cluster_rows = T,
              heatmap_legend_param = list(at = c(0,1), labels = c("Missing value", "Valid value")))
dev.off()

## PD1 exp for each antibody condition ####
pdf("~/ip_proteomics/figures/pd1exp_condition.pdf", width = 7, height = 5)
ggplot(annotation_df, aes(x=antibody, y=pd1_exp)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=2) + theme(text = element_text(size = 8))
dev.off()

## visualize number of proteins/sample ####
protein_data_number <- ifelse(is.na(protein_data), 0, 1)
protein_data_number <- as.data.frame(colSums(protein_data_number))
protein_data_number$experiment <- rownames(protein_data_number)
keep <- protein_data_number[apply(protein_data_number, MARGIN = 1, function(x) !all(is.na(x) == TRUE)), ]
protein_data_number <- cbind(protein_data_number, annotation_df[,c(2,3,4), drop =FALSE])
protein_data_number <- protein_data_number %>% arrange(desc(abs(`colSums(protein_data_number)`)))


pdf("/home/degan/ip_proteomics/figures/proteins_persample.pdf", width = 5, height = 5)
ggplot(protein_data_number, aes(x = experiment, y = `colSums(protein_data_number)`, fill = antibody)) + geom_bar(stat = "identity") + 
  geom_hline(yintercept = nrow(keep)) + labs(title = "Proteins per sample", 
                                             x = "", y = "Number of proteins")  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4)) 
dev.off()

## Filter proteins ####
## Filter for proteins that are identified in 8 out of 8 replicates of at least one condition ####
filter_missval_own <- function(protein.df, exp_design.df, thr = 0) {
  
  # Make assay values binary (1 = valid value)
  bin_data <- protein.df
  idx <- is.na(protein.df)
  bin_data[!idx] <- 1
  bin_data[idx] <- 0
  
  # Filter se on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(ID, value, -rowname) %>%
    left_join(., data.frame(exp_design.df), by = c("ID" = "label")) %>%
    group_by(rowname, condition) %>%
    summarize(miss_val = n() - sum(value)) %>%
    filter(miss_val <= thr) %>%
    spread(condition, miss_val)
  protein_fltrd <- protein.df[keep$rowname, ]
  return(protein_fltrd)
}

protein_data_filtered <- filter_missval_own(protein_data, experimental_design, thr=0)

## Plot a barplot of the protein identification overlap between samples ####

# Make a binary long data.frame (1 = valid value, 0 = missing value)
protein_data_number <- protein_data_filtered %>%
  rownames_to_column() %>%
  gather(ID, bin, -rowname) %>%
  mutate(bin = ifelse(is.na(bin), 0, 1))
# Identify the number of experiments a protein was observed
stat <- protein_data_number %>%
  group_by(rowname) %>%
  summarize(sum = sum(bin))
# Get the frequency of the number of experiments proteins
# were observerd and plot these numbers
table <- table(stat$sum) %>% data.frame()
pdf("/home/degan/ip_proteomics/figures/protein_sample_overlap.pdf", height = 5, width = 7)
ggplot(table, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_col() +
  scale_fill_grey(start = 0.8, end = 0.2) +
  labs(title = "Protein identifications overlap",
       x = "Identified in number of samples",
       y = "Number of proteins") +
  theme(legend.position="none")
dev.off()

## normalize filtered protein data using global vsn transformation ####
fit <- vsn::vsn2(as.matrix(protein_data_filtered))
meanSdPlot(fit) + theme_bw()
protein_data_norm <- predict(fit, newdata= as.matrix(protein_data_filtered))

## qq plots for before and after normalization ####
###############################################################################

qqnorm(protein_data_filtered[1,], pch = 1, frame = FALSE)
qqline(protein_data_filtered[1,], col = "steelblue", lwd = 2)

qqnorm(protein_data_norm[1,], pch = 1, frame = FALSE)
qqline(protein_data_norm[1,], col = "steelblue", lwd = 2)

## Horizontal Imputation ####
# GroupMeanGD
# new version
# group mean imputation with normal distribution correction
# or stochastic group mean value imputation
HorizontalInputation <- function (df) 
{
  nSamples = dim(df)[2]
  nFeatures = dim(df)[1]
  count.NAs = apply(!is.na(df), 1, sum)
  dataSet.filtered = df[which(count.NAs >= 4), ]
  dataSet.not_filtered = df[which(count.NAs < 4), ]
  mean.samples = apply(dataSet.filtered, 1, mean, na.rm = T)
  protSD = apply(dataSet.filtered, 1, sd, na.rm = T)
  for (i in 1:(nrow(dataSet.filtered))) {
    dataSet.filtered.temp = rnorm(nSamples, mean = mean.samples[i], sd = protSD[i])
    dataSet.filtered[i,which(is.na(dataSet.filtered[i,]))] = dataSet.filtered.temp[which(is.na(dataSet.filtered[i,]))]
  }
  #k <- which(is.na(dataSet.imputed), arr.ind=TRUE)
  #dataSet.imputed[k] <- rowMeans(dataSet.imputed, na.rm=TRUE)[k[,1]]
  dataSet.final <- rbind(dataSet.filtered,dataSet.not_filtered)
  dataSet.final <- dataSet.final[order(row.names(dataSet.final)), ]
  return(dataSet.final)
}

## calling function ###
groupsP <- levels(as.factor(experimental_design$condition))
all <- HorizontalInputation(protein_data_norm[,which(experimental_design$condition %in% groupsP[1])])
for (i in 2:length(groupsP)) {
  test1 <- HorizontalInputation(protein_data_norm[,which(experimental_design$condition %in% groupsP[i])])
  all <- cbind(all,test1)
}

## Horizontal imputation + MinProb ####
impdata.Horizontal_MinProb <- impute.MinProb(as.matrix(all), q = 0.01, tune.sigma = 1)

## save imputed protein data matrix ####
saveRDS(impdata.Horizontal_MinProb, "/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds")
saveRDS(annotation_df, "/home/degan/ip_proteomics/inputs/annotation_df.Rds")

write.csv(data.frame(impdata.Horizontal_MinProb), "/home/degan/ip_proteomics/inputs/imputed_protein_matrix.csv", row.names = T)
write.csv(annotation_df, "/home/degan/ip_proteomics/inputs/annotation_df.csv", row.names = T)

## visualize number of proteins/sample ####
protein_data_number <- ifelse(is.na(impdata.Horizontal_MinProb), 0, 1)
protein_data_number <- as.data.frame(colSums(protein_data_number))
protein_data_number$experiment <- rownames(protein_data_number)
keep <- protein_data_number[apply(protein_data_number, MARGIN = 1, function(x) !all(is.na(x) == TRUE)), ]
protein_data_number <- cbind(protein_data_number, annotation_df[,c(2,3,4), drop =FALSE])
protein_data_number <- protein_data_number %>% arrange(desc(abs(`colSums(protein_data_number)`)))

ggplot(protein_data_number, aes(x = experiment, y = `colSums(protein_data_number)`, fill = antibody)) + geom_bar(stat = "identity") + 
  geom_hline(yintercept = nrow(keep)) + labs(title = "Proteins per sample", 
                                             x = "", y = "Number of proteins")  + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4)) 



