#Heatmap - Fig. 1A
library("readr")
load("CD4_RGSS_data_aftercombatM_case.RData")
clinical = read.csv("CD4_RGSS_clinical_case.csv", header = T)
Group_Disease_2 = clinical$Disease_Type
pvalue =apply(CD4_RGSS_data_aftercombatM_case, 1, function(y){return(summary(aov(y~as.factor(Group_Disease_2)))[[1]][5][1,1])})
pvalue.fdr = p.adjust(pvalue,method = "fdr")
pvalue.bonferroni = p.adjust(pvalue,method = "bonferroni")
cg.Group_Disease = data.frame(pvalue,pvalue.fdr,pvalue.bonferroni)
cg.Group_Disease = cg.Group_Disease[order(cg.Group_Disease[,1]),]
anno = read_tsv("HM450KAnnotation.txt")
seq = match(rownames(cg.Group_Disease),anno$ID);
anno_top = anno[seq,c(11,12,21,23,25)];
cg.Group_Disease = data.frame(cg.Group_Disease, anno_top);
write.csv(cg.Group_Disease,file="cg.Group_Disease_ANOVA_case_aftercombatM.csv",quote=F)
ANOVA = as.data.frame(read_csv("cg.Group_Disease_ANOVA_case_aftercombatM.csv", col_names = T))
sig.cgs = ANOVA[1:50,1]
cgs.idx = match(sig.cgs, rownames(CD4_RGSS_data_aftercombatM_case))
data.heatmap = data.matrix(CD4_RGSS_data_aftercombatM_case[cgs.idx,])
sample.info = read.csv("CD4_RGSS_clinical_case.csv", header = T)
Disease = sample.info$Disease_Type
library(gplots)
library(RColorBrewer)
col = rev(colorRampPalette(brewer.pal(10,'RdBu'))(250))
sample.color = ifelse(Disease =="GD","red", ifelse(Disease=="RA","violet", ifelse(Disease=="SLE","cornflowerblue","cyan")))
heatmap.2(data.heatmap, col=col, scale="row", ColSideColors=sample.color, breaks = seq(-5,5,0.04), Rowv = F, dendrogram = "column", 
          key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none",labRow = F, labCol = F)
legend("topleft",legend = c("GD","RA","SLE","SSc"), fill = c("red","violet","cornflowerblue","cyan"), bty = "n", cex = 0.7)
#Heatmap - Fig. 1B
library("readr")
load("CD4_RGSS_data_aftercombatM.RData")
comparison = as.data.frame(read_csv("cg.comparison_aftercombatM.csv", col_names = T))
sig.cgs = comparison[1:50,1]
cgs.idx = match(sig.cgs, rownames(CD4_RGSS_data_aftercombatM))
sample.info = read.csv("CD4_RGSS_clinical_new2.csv", header = T)
data.heatmap = data.matrix(CD4_RGSS_data_aftercombatM[cgs.idx,])
Group_Disease_2 = sample.info$Group_Disease_2
library(gplots)
library(RColorBrewer)
col = rev(colorRampPalette(brewer.pal(10,'RdBu'))(250))
sample.color = ifelse(Group_Disease_2 == "Control","blue",ifelse(Group_Disease_2 =="GD","red", ifelse(Group_Disease_2=="RA","violet", 
                                                                                                      ifelse(Group_Disease_2=="SLE","cornflowerblue","cyan"))))
heatmap.2(data.heatmap, col=col, scale="row", ColSideColors=sample.color, breaks = seq(-5,5,0.04), Rowv = F, dendrogram = "column", 
          key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none",labRow = F, labCol = F)
legend("topleft",legend = c("GD","RA","SLE","SSc","Control"), fill = c("red","violet","cornflowerblue","cyan","blue"), 
       bty = "n", cex = 0.7)
#Heatmap - Supplementary Fig. 3
library("readr")
load("CD4_RGSS_data_aftercombatM.RData")
ANOVA = as.data.frame(read_csv("cg.Group_Disease_ANOVA_5groups_aftercombatM.csv", col_names = T))
sig.cgs = ANOVA[1:50,1]
cgs.idx = match(sig.cgs, rownames(CD4_RGSS_data_aftercombatM))
sample.info = read.csv("CD4_RGSS_clinical_new2.csv", header = T)
data.heatmap = data.matrix(CD4_RGSS_data_aftercombatM[cgs.idx,])
Group_Disease_2 = sample.info$Group_Disease_2
library(gplots)
library(RColorBrewer)
col = rev(colorRampPalette(brewer.pal(10,'RdBu'))(250))
sample.color = ifelse(Group_Disease_2 == "Control","blue",ifelse(Group_Disease_2 =="GD","red", ifelse(Group_Disease_2=="RA","violet", 
                      ifelse(Group_Disease_2=="SLE","cornflowerblue","cyan"))))
heatmap.2(data.heatmap, col=col, scale="row", ColSideColors=sample.color, breaks = seq(-5,5,0.04), Rowv = F, dendrogram = "column", 
          key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none",labRow = F, labCol = F)
legend("topleft",legend = c("GD","RA","SLE","SSc","Control"), fill = c("red","violet","cornflowerblue","cyan","blue"), 
       bty = "n", cex = 0.7)
#Heatmap - Fig. 3A
load("CD4_RGSS_data_aftercombatM.RData")
comparison = read.csv("cg.comparison_aftercombatM.csv", header = T,stringsAsFactors = F)
DMS = comparison[which(comparison[,3] < 0.01),]
seq1 = which(DMS$UCSC_RefGene_Name == "MX1")
MX1 = DMS[seq1,]
MX1_beta = as.data.frame(CD4_RGSS_data_aftercombatM[MX1$X,])
for(i in 1:length(rownames(MX1_beta))){
  rownames(MX1_beta)[i] = paste(rownames(MX1_beta)[i],"(MX1)",sep = " ")}
seq2 = which(DMS$UCSC_RefGene_Name == "OAS1")
OAS1 = DMS[seq2,]
OAS1_beta = as.data.frame(CD4_RGSS_data_aftercombatM[OAS1$X,])
for(i in 1:length(rownames(OAS1_beta))){
  rownames(OAS1_beta)[i] = paste(rownames(OAS1_beta)[i],"(OAS1)",sep = " ")}
seq3 = which(DMS$UCSC_RefGene_Name == "USP18")
USP18 = DMS[seq3,]
USP18_beta = as.data.frame(CD4_RGSS_data_aftercombatM[USP18$X,])
for(i in 1:length(rownames(USP18_beta))){
  rownames(USP18_beta)[i] = paste(rownames(USP18_beta)[i],"(USP18)",sep = " ")}
seq4 = which(DMS$UCSC_RefGene_Name == "IFIT1")
IFIT1 = DMS[seq4,]
IFIT1_beta = as.data.frame(CD4_RGSS_data_aftercombatM[IFIT1$X,])
for(i in 1:length(rownames(IFIT1_beta))){
  rownames(IFIT1_beta)[i] = paste(rownames(IFIT1_beta)[i],"(IFIT1)",sep = " ")}
seq5 = which(DMS$UCSC_RefGene_Name == "IRF7")
IRF7 = DMS[seq5,]
IRF7_beta = as.data.frame(CD4_RGSS_data_aftercombatM[IRF7$X,])
for(i in 1:length(rownames(IRF7_beta))){
  rownames(IRF7_beta)[i] = paste(rownames(IRF7_beta)[i],"(IRF7)",sep = " ")}
seq6 = which(DMS$UCSC_RefGene_Name == "RSAD2")
RSAD2 = DMS[seq6,]
RSAD2_beta = as.data.frame(CD4_RGSS_data_aftercombatM[RSAD2$X,])
for(i in 1:length(rownames(RSAD2_beta))){
  rownames(RSAD2_beta)[i] = paste(rownames(RSAD2_beta)[i],"(RSAD2)",sep = " ")}
data.heatmap = data.matrix(rbind(MX1_beta,OAS1_beta,USP18_beta,IFIT1_beta,IRF7_beta,RSAD2_beta))
clinical = read.csv("CD4_RGSS_clinical.csv", header = T)
Sample_Group = clinical$Sample_Group
sample.color = ifelse(Sample_Group == "Case","indianred1","cornflowerblue")
library(gplots)
library(RColorBrewer)
col = rev(colorRampPalette(brewer.pal(10,'RdBu'))(250))
heatmap.2(data.heatmap, col=col, scale="row", ColSideColors=sample.color, breaks = seq(-5,5,0.04), Rowv = F, dendrogram = "column", 
          key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none",labCol = F, margins = c(5,10))
legend("topleft",legend = c("GD/RA/SLE/SSc Patients","Control"), fill = c("indianred1","cornflowerblue"), bty = "n", cex = 0.7)