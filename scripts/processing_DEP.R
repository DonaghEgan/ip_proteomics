## Library ####
library(DEP)

## action points ####

# remove lowly abundant proteins: must be present in all 3 replicates of atleast 1 condition
# compare pull downs to whole proteome - check if significant proteins are correlated with
# most abundant proteins. 


## Reading in file ####
pgroups <- read.table(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/proteinGroups_PPI_Martina_Jan.txt", 
                      header = T, sep = "\t",quote='',stringsAsFactors = FALSE,comment.char="")

## Removing contaminants and unnecessary characters ####
pgroups <- pgroups[pgroups$Reverse != "+" & pgroups$Potential.contaminant != "+", ]

## Make unique names using the annotation in the "Gene.names" column as primary ####
## names and the annotation in "Protein.IDs" as name for those that do not have an gene name.####
pgroups_unique <- make_unique(pgroups, "Gene.names", "Protein.IDs", delim = ";")
pgroups_unique$name %>% duplicated() %>% any()

## Generate SummarizedExperiment Object ####
LFQ_columns <- grep("LFQ", colnames(pgroups_unique)) # get LFQ column numbers

assay_p <- data.frame(pgroups_unique[,grep("LFQ", colnames(pgroups_unique))], row.names = pgroups_unique$name)
rowdata_p <- data.frame(row.names = pgroups_unique$name, pgroups_unique$name) 

col_data_p <- data.frame(condition = ifelse(grepl("Con",colnames(protein_data)), "Control", "Antibody"), 
                         row.names =  colnames(pgroups_unique),
                         antibody =  sapply(strsplit(colnames(protein_data), "_"), "[", 1),
                         time = sapply(strsplit(colnames(protein_data), "_"), "[", 2),
                         batch = sapply(strsplit(colnames(protein_data), "_"), "[", 4),
                         pd1_exp = as.numeric(protein_data["Q15116",]))

se_proteins <- SummarizedExperiment::SummarizedExperiment(assays = assay_p, rowData = rowdata_p,
                                                          colData = col_data_p)

## filter on missing values ####

## Plot a barplot of the protein identification overlap between samples ####
plot_frequency(se_proteins)

?makeSummarizedExperimentFromDataFrame

## Filter for proteins that are identified in 2 out of 3 replicates of at least one condition ####
data_filter <- filter_missval(se_proteins, thr = 1)
plot_numbers(data_filter)

