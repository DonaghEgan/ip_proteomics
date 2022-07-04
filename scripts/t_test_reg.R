# t-test from regression coefficients yo!

ttest_reg <- function(y1, x1, x2, y2, x3, y3){
  
  # perform regression 1
  reg1 = glm(y1~x1)
  # perform regression 2. 
  reg2 = glm(y2~x2)
  # perform regression 3. 
  reg3 = glm(y3~x3)
  
  # get t statistic 
  #[ formula t = (b1 - b2)/(se(b))]
  
  # get coefficients. 
  b1 = reg1$coefficients[2]
  b2 = reg2$coefficients[2]
  b3 = reg3$coefficients[2]
  
  # pool variance. 
  se_b1 = summary(reg1)$coefficients[2,2]
  se_b2 = summary(reg2)$coefficients[2,2]
  se_b3 = summary(reg3)$coefficients[2,2]
  
  # get t-statistic
  n1 = length(y1)
  n2 = length(y2)
  n3 = length(y3)
  
  P_val = ttest_two(b1, b2, se_b1, se_b2, n1, n2)
  return(P_val)}
  
ttest_two <- function(b1, b2, se_b1, se_b2, n1 ,n2){
  t = (b1 -b2)/sqrt((se_b1^2/n1) + (se_b2^2/n2))
  t = as.numeric(t)
  deg_free = n1+n2 -2
  # calculate t-test
  P_val_1 = pt(q= t, df = deg_free) # b1 less than b2
  P_val_2 = 1-pt(t, deg_free) #  b1 greater than b2
  P_Val = 2*min(P_val_2, P_val_1)
  return(P_Val)
  
}

test_plot <- function(b1, b2, se_b1, se_b2)
  
  coef_df <- melt(data.frame(csa = b1, niv = b2, pem=b3))
  names(coef_df)[1] <- "antibody"
  names(coef_df)[2] <- "coef_val"
  
  se_df <- melt(data.frame(csa = se_b1, niv = se_b2, pem=se_b3))
  names(se_df)[1] <- "antibody"
  names(se_df)[2] <- "se_val"
  
  res_df <- merge(coef_df, se_df)
  
  pdf("/home/degan/ip_proteomics/figures/UMAP/coef_antibody_mcm_pd1.pdf", height = 3, width = 3)
  ggplot(res_df, aes(x=antibody, y=coef_val, ymin = coef_val-se_val, ymax = coef_val+se_val, color=antibody)) + 
  geom_errorbar(width = 0.2) + geom_point(size = 1.5) + theme_classic() + ylab("Coef") + xlab("")  + ylim(0,1)
  dev.off()











