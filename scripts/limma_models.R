## Library ####
library(limma)
library(data.table)

## read in exp design - created by Vadim ####
experimental_design <- read.csv(file = "/mnt/Data/Vadim/POIAZ/Vadim/PPI/PPI_Martina_Jan/experimental_design1.csv", sep = "\t", 
                                colClasses=c("character","character"), header = F)
colnames(experimental_design) <- c("label","condition")

##remove stimulation from exp design
experimental_design$stimulation <- ifelse(grepl("Stim",experimental_design$label), "Stimulated", "NotStim")
experimental_design <- experimental_design[which(experimental_design$stimulation == "NotStim"),]


## creating design matrix ####
groups <- as.factor(experimental_design$condition)
designlimmaP <- as.data.frame(model.matrix(~ 0 + groups))
names(designlimmaP) <- gsub("groups", "", names(designlimmaP))

## contrast each control and PD1 for each antibody and timepoint  ####  
contrastP <- limma::makeContrasts(paste(levels(groups)[16],"-",levels(groups)[1],sep = ""),
                    paste(levels(groups)[17],"-",levels(groups)[2],sep = ""),
                    paste(levels(groups)[18],"-",levels(groups)[3],sep = ""),
                    paste(levels(groups)[19],"-",levels(groups)[4],sep = ""),
                    paste(levels(groups)[20],"-",levels(groups)[5],sep = ""),
                    
                    
                    paste(levels(groups)[21],"-",levels(groups)[6],sep = ""),
                    paste(levels(groups)[22],"-",levels(groups)[7],sep = ""),
                    paste(levels(groups)[23],"-",levels(groups)[8],sep = ""),
                    paste(levels(groups)[24],"-",levels(groups)[9],sep = ""),
                    paste(levels(groups)[25],"-",levels(groups)[10],sep = ""),
                    
                    paste(levels(groups)[26],"-",levels(groups)[11],sep = ""),
                    paste(levels(groups)[27],"-",levels(groups)[12],sep = ""),
                    paste(levels(groups)[28],"-",levels(groups)[13],sep = ""),
                    paste(levels(groups)[29],"-",levels(groups)[14],sep = ""),
                    paste(levels(groups)[30],"-",levels(groups)[15],sep = ""),
                    levels=designlimmaP)
contrastP <- t(na.omit(t(contrastP)))

fit_ctrl_ant <- lmFit(imputed_nostim, designlimmaP) # to check
fit_ctrl_ant <- contrasts.fit(fit_ctrl_ant, contrastP)
fit_ctrl_ant<- eBayes(fit_ctrl_ant) 

coef_ctrl_ant <- list() ### varian of code: https://www.r-bloggers.com/concatenating-a-list-of-data-frames/
for (s in 1:length(colnames(contrastP))) {
  coef_ctrl_ant[[s]] <- topTable(fit_ctrl_ant, coef=s, n=Inf)
  setDT(coef_ctrl_ant[[s]], keep.rownames = TRUE)[]
  coef_ctrl_ant[[s]]$type <- colnames(contrastP)[s]
}

diff_exp_protein <- do.call(rbind, coef_ctrl_ant)




