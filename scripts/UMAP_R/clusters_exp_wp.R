library(RColorBrewer)

## ip proteomics ####
cluster_info <- readRDS("/home/degan/ip_proteomics/inputs/cluster_assignments.Rds")
annotation_df <- readRDS("/home/degan/ip_proteomics/inputs/annotation_df.Rds")
imputed_data <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/imputed_protein_matrix.Rds"))
gene_names <- data.frame(readRDS("/home/degan/ip_proteomics/inputs/pggroups.meta.Rds"))

## whole proteomics ####
## Read in imputed data matrix and anntotations ####
w_proteomics <- data.frame(readRDS("/home/degan/whole_proteomics/saved_objects/processed_proteomics.Rds"))
exp_design_wp <- readRDS("/home/degan/whole_proteomics/saved_objects/experimental_design.Rds")

## Phospho proteomics ####
## Read in imputed data matrix and anntotations ####
phos_proteomics <- readRDS("/mnt/Data/Vadim/POIAZ/Donagh/proteomics22_unified/processing_&_normalisation/saved_objects/phos_imputed_for_int.Rds")
exp_design_wp <- readRDS("/home/degan/whole_proteomics/saved_objects/experimental_design.Rds")


## CALCULATE AVERAGE DIFFERENCE BETWEEN PD1 AND CTRL FOR EACH CLUSTERS - WP ####
wp_diff <- list()
for (i in unique(cluster_info$cluster_id)) {
  ## proteins for given cluster
  cluster_proteins <- cluster_info[cluster_info$cluster_id == i,]
  ## remove semi colon from proteins - not present in whole proteomics
  cluster_proteins <- strsplit(cluster_proteins$Row.names, ";", fixed=F)
  cluster_proteins <- Reduce(c,cluster_proteins)
  
  ## subset wp based on cluster proteins
  cluster_wp <- w_proteomics[rownames(w_proteomics) %in% cluster_proteins,]
  cluster_wp <- t(cluster_wp)
  
  ## add sample info 
  cluster_wp <- merge(cluster_wp, exp_design_wp, by = 0)
  cluster_wp <- cluster_wp %>% column_to_rownames("Row.names")
  
  ## mean proteins abundance for pd1 and ctrl 
  cluster_wp_mean <- aggregate(.~ type, data = cluster_wp[,c(1:(ncol(cluster_wp)-7),(ncol(cluster_wp) - 2))], mean)
  cluster_wp_mean <- cluster_wp_mean %>% column_to_rownames("type") %>% t() %>% data.frame()
  
  ## calculate diff 
  mean_diff <- mean(cluster_wp_mean$PD1) - mean(cluster_wp_mean$CNTR)
  
  ## add to list 
  wp_diff[[i]] <- mean_diff 
}

## CALCULATE AVERAGE DIFFERENCE BETWEEN PD1 AND CTRL FOR EACH CLUSTERS - IP ####

ip_diff <- list()
for (i in unique(cluster_info$cluster_id)) {
  
  ## proteins for given cluster
  cluster_proteins <- cluster_info[cluster_info$cluster_id == i,]

  ## subset ip based on cluster proteins
  cluster_ip <- imputed_data[rownames(imputed_data) %in% cluster_proteins$Row.names,]
  cluster_ip <- t(cluster_ip)
  
  ## add sample info 
  cluster_ip <- merge(cluster_ip, annotation_df, by = 0)
  cluster_ip <- cluster_ip %>% column_to_rownames("Row.names")
  
  ## mean proteins abundance for pd1 and ctrl 
  cluster_ip_mean <- aggregate(.~ condition, data = cluster_ip[,c(1:(ncol(cluster_ip)-5))], mean)
  cluster_ip_mean <- cluster_ip_mean %>% column_to_rownames("condition") %>% t() %>% data.frame()
  
  ## calculate diff 
  mean_diff <- mean(cluster_ip_mean$PD1) - mean(cluster_ip_mean$Con)
  
  ## add to list 
  ip_diff[[i]] <- mean_diff 
}

## COMBINE to dataframe for comparison
ip_wp_df <- data.frame(row.names = names(ip_diff), ip = unlist(ip_diff), wp = unlist(wp_diff))
ip_wp_df$complex <- rownames(ip_wp_df)

ggplot(ip_wp_df, aes(ip, wp)) + geom_point() + geom_smooth(method = "lm")

## PLOT RELATIONSHIP ####
################################################################################

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pdf("/home/degan/ip_proteomics/figures/UMAP/relationship_diff_complex_wp_ip.pdf", width = 4, height = 4)
ggscatter(ip_wp_df, x = "wp", y = "ip", color = "complex",
          add = "reg.line", conf.int = FALSE, size= 0.9, label.rectangle = T,
          cor.coef = TRUE, cor.method = "pearson", add.params = list(color = pal_npg("nrc")(10)[1], linetype = "dashed", size = 0.6, alpha = 0.7),
          xlab = "IP: Mean Difference PD1 vs Ctr ", ylab = "WP: Mean Difference PD1 vs Ctr") +
          theme(axis.title = element_text(size = 9),
          legend.position = "bottom", legend.text = element_text(size=8)) + 
          scale_color_manual(values =  col_vector)
dev.off()

## CALCULATE AVERAGE DIFFERENCE BETWEEN PD1 AND CTRL FOR EACH CLUSTERS - PHOSPHO ####
phos_diff <- list()
for (i in unique(cluster_info$cluster_id)) {
  ## proteins for given cluster
  cluster_proteins <- cluster_info[cluster_info$cluster_id == i,]
  ## remove semi colon from proteins - not present in whole proteomics
  cluster_proteins <- strsplit(cluster_proteins$Row.names, ";", fixed=F)
  cluster_proteins <- Reduce(c,cluster_proteins)
  
  ## subset wp based on cluster proteins
  cluster_phos <- w_proteomics[rownames(w_proteomics) %in% cluster_proteins,]
  cluster_wp <- t(cluster_wp)
  
  ## add sample info 
  cluster_wp <- merge(cluster_wp, exp_design_wp, by = 0)
  cluster_wp <- cluster_wp %>% column_to_rownames("Row.names")
  
  ## mean proteins abundance for pd1 and ctrl 
  cluster_wp_mean <- aggregate(.~ type, data = cluster_wp[,c(1:(ncol(cluster_wp)-7),(ncol(cluster_wp) - 2))], mean)
  cluster_wp_mean <- cluster_wp_mean %>% column_to_rownames("type") %>% t() %>% data.frame()
  
  ## calculate diff 
  mean_diff <- mean(cluster_wp_mean$PD1) - mean(cluster_wp_mean$CNTR)
  
  ## add to list 
  wp_diff[[i]] <- mean_diff 
}

