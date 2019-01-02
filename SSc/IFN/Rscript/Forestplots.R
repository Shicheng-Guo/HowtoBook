#Forestplot - Fig. 2
library(readr)
load("CD4_RGSS_data_beforecombat.RData")
CD4_DMS = read_csv("cg.comparison_aftercombatM.csv")
CD4_DMS_50 = CD4_DMS[1:50,]
seq = match(CD4_DMS_50$X1,rownames(CD4_RGSS_data_beforecombat))
CD4_DMS50_beta = CD4_RGSS_data_beforecombat[seq,]
CD4_clinical = read_csv("CD4_RGSS_clinical.csv")
case="Case"
control="Control"
Diseases = unique(CD4_clinical$Disease_Type)
for(i in 1:length(Diseases)){
  seq_case = ifelse(CD4_clinical$Sample_Group==case,ifelse(CD4_clinical$Disease_Type==Diseases[i],TRUE,FALSE),FALSE)
  seq_control = ifelse(CD4_clinical$Sample_Group==control,ifelse(CD4_clinical$Disease_Type==Diseases[i],TRUE,FALSE),FALSE)
  total_case = apply(CD4_DMS50_beta,1,function(x){return(length(x[seq_case]))})
  mean_case = apply(CD4_DMS50_beta,1,function(x){return(mean(x[seq_case],na.rm = T))})
  sd_case = apply(CD4_DMS50_beta,1,function(x){return(sd(x[seq_case],na.rm = T))})
  total_control = apply(CD4_DMS50_beta,1,function(x){return(length(x[seq_control]))})
  mean_control = apply(CD4_DMS50_beta,1,function(x){return(mean(x[seq_control],na.rm = T))})
  sd_control = apply(CD4_DMS50_beta,1,function(x){return(sd(x[seq_control],na.rm = T))})
  table = data.frame(total_case,mean_case,sd_case,total_control,mean_control,sd_control)
  write.csv(table, file = paste("CD4_DMSbeta_", Diseases[i], ".csv",sep=""), quote=F, row.names = T)
}
CD4_GD = read_csv("CD4_DMSbeta_GD.csv")
CD4_RA = read_csv("CD4_DMSbeta_RA.csv")
CD4_SLE = read_csv("CD4_DMSbeta_SLE.csv")
CD4_SSc = read_csv("CD4_DMSbeta_SSc.csv")
library(metafor)
par(mfrow = c(3,2))
CD4_DMS_9 = rbind(CD4_GD[9,],CD4_RA[9,],CD4_SLE[9,],CD4_SSc[9,])
CD4_DMS_9$X1 = c("GD","RA","SLE","SSc")
colnames(CD4_DMS_9)[1] = "Disease"
meta_CD4_DMS9 = rma.uni(n1i = total_case, n2i = total_control, m1i = mean_case, m2i = mean_control, sd1i = sd_case, sd2i = sd_control, 
                        data=CD4_DMS_9, measure = "SMD", method="DL")
forest_CD4_DMS9 = forest(meta_CD4_DMS9,slab = paste(CD4_DMS_9$Disease," vs. control", sep = ""), mlab="Random-effect Model")
CD4_DMS_20 = rbind(CD4_GD[20,],CD4_RA[20,],CD4_SLE[20,],CD4_SSc[20,])
CD4_DMS_20$X1 = c("GD","RA","SLE","SSc")
colnames(CD4_DMS_20)[1] = "Disease"
meta_CD4_DMS20 = rma.uni(n1i = total_case, n2i = total_control, m1i = mean_case, m2i = mean_control, sd1i = sd_case, sd2i = sd_control, 
                         data=CD4_DMS_20, measure = "SMD", method="DL")
forest_CD4_DMS20 = forest(meta_CD4_DMS20,slab = paste(CD4_DMS_20$Disease," vs. control", sep = ""), mlab="Random-effect Model")
CD4_DMS_11 = rbind(CD4_GD[11,],CD4_RA[11,],CD4_SLE[11,],CD4_SSc[11,])
CD4_DMS_11$X1 = c("GD","RA","SLE","SSc")
colnames(CD4_DMS_11)[1] = "Disease"
meta_CD4_DMS11 = rma.uni(n1i = total_case, n2i = total_control, m1i = mean_case, m2i = mean_control, sd1i = sd_case, sd2i = sd_control, 
                         data=CD4_DMS_11, measure = "SMD", method="DL")
forest_CD4_DMS11 = forest(meta_CD4_DMS11,slab = paste(CD4_DMS_11$Disease," vs. control", sep = ""), mlab="Random-effect Model")
CD4_DMS_17 = rbind(CD4_GD[17,],CD4_RA[17,],CD4_SLE[17,],CD4_SSc[17,])
CD4_DMS_17$X1 = c("GD","RA","SLE","SSc")
colnames(CD4_DMS_17)[1] = "Disease"
meta_CD4_DMS17 = rma.uni(n1i = total_case, n2i = total_control, m1i = mean_case, m2i = mean_control, sd1i = sd_case, sd2i = sd_control, 
                         data=CD4_DMS_17, measure = "SMD", method="DL")
forest_CD4_DMS17 = forest(meta_CD4_DMS17,slab = paste(CD4_DMS_17$Disease," vs. control", sep = ""), mlab="Random-effect Model")
CD4_DMS_16 = rbind(CD4_GD[16,],CD4_RA[16,],CD4_SLE[16,],CD4_SSc[16,])
CD4_DMS_16$X1 = c("GD","RA","SLE","SSc")
colnames(CD4_DMS_16)[1] = "Disease"
meta_CD4_DMS16 = rma.uni(n1i = total_case, n2i = total_control, m1i = mean_case, m2i = mean_control, sd1i = sd_case, sd2i = sd_control, 
                         data=CD4_DMS_16, measure = "SMD", method="DL")
forest_CD4_DMS16 = forest(meta_CD4_DMS16,slab = paste(CD4_DMS_16$Disease," vs. control", sep = ""), mlab="Random-effect Model")
CD4_DMS_33 = rbind(CD4_GD[33,],CD4_RA[33,],CD4_SLE[33,],CD4_SSc[33,])
CD4_DMS_33$X1 = c("GD","RA","SLE","SSc")
colnames(CD4_DMS_33)[1] = "Disease"
meta_CD4_DMS33 = rma.uni(n1i = total_case, n2i = total_control, m1i = mean_case, m2i = mean_control, sd1i = sd_case, sd2i = sd_control, 
                         data=CD4_DMS_33, measure = "SMD", method="DL")
forest_CD4_DMS33 = forest(meta_CD4_DMS33,slab = paste(CD4_DMS_33$Disease," vs. control", sep = ""), mlab="Random-effect Model")
#Forestplot - Supplementary Fig. 4
CD4_DMS_42 = rbind(CD4_GD[42,],CD4_RA[42,],CD4_SLE[42,],CD4_SSc[42,])
CD4_DMS_42$X1 = c("GD","RA","SLE","SSc")
colnames(CD4_DMS_42)[1] = "Disease"
meta_CD4_DMS42 = rma.uni(n1i = total_case, n2i = total_control, m1i = mean_case, m2i = mean_control, sd1i = sd_case, sd2i = sd_control, 
                         data=CD4_DMS_42, measure = "SMD", method="DL")
forest_CD4_DMS42 = forest(meta_CD4_DMS42,slab = paste(CD4_DMS_42$Disease," vs. control", sep = ""), mlab="Random-effect Model")