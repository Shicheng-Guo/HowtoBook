#Barplot - Supplementary Fig. 1

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ggsci)
library(ggthemes)
##CGI feature
#DMS
CD4 = read.csv("cg.comparison_aftercombatM.csv", header = T, stringsAsFactors = F)
CD4_DMS = CD4[which(CD4[,3] < 0.01),]
CGI = CD4_DMS$Relation_to_UCSC_CpG_Island
CGI[is.na(CGI)] = "Open_Sea"
CGI = as.data.frame(CGI)
CGIprop = as.data.frame(prop.table(table(CGI)))
CGIprop$Methylation_Sites = "Differential"
#All data
load("CD4_RGSS_data_aftercombatM.RData")
library(readr)
anno = read_tsv("HM450KAnnotation.txt")
seq = match(rownames(CD4_RGSS_data_aftercombatM),anno$ID);
anno_all = anno[seq,c(11,12,21,23,25)];
anno_all = as.data.frame(anno_all)
rownames(anno_all) = row.names(CD4_RGSS_data_aftercombatM)
write.csv(anno_all,file="anno_aftercombatM.csv",quote=F)
CGI_all = anno_all$Relation_to_UCSC_CpG_Island
CGI_all[is.na(CGI_all)] = "Open_Sea"
CGI_all = as.data.frame(CGI_all)
write.csv(CGI_all,file="CD4_CpG_subregion_all.csv",quote=F)
CGIprop_all = as.data.frame(prop.table(table(CGI_all)))
CGIprop_all$Methylation_Sites = "All"
colnames(CGIprop_all)[1] = "CGI"
CGIprop_merge = rbind(CGIprop,CGIprop_all)
colnames(CGIprop_merge)[2] = "Proportion"
barplot = ggplot(CGIprop_merge, aes(x = CGI, y = Proportion, fill = Methylation_Sites)) + geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 0.6, 0.1),limits = c(0,0.6)) +
  theme_classic()+
  theme(legend.position=c(0.88,0.90))+theme(legend.text = element_text(colour="black", size = 11))+theme(legend.background = element_rect(colour = "black"))+
  theme(axis.text.x=element_text(colour = "black", size=11))+ theme(axis.text.y=element_text(colour = "black", size=11))+
  theme(axis.title.x=element_blank()) +  theme(axis.title.y= element_text(colour = "black", size=12))+
  scale_color_npg()+scale_fill_npg()
ggsave(filename = "CGIfeature_CD4_dva_aftercombatM_001_npg.tiff", plot = barplot, width = 7, height = 6, dpi = 500)
ggsave(filename = "CGIfeature_CD4_dva_aftercombatM_001_npg.pdf", plot = barplot, width = 7, height = 6, dpi = 500)

##Gene feature
#DMS
Summary = read.csv("Genefeature_DMS_aftercombatM_001.csv", header = T, stringsAsFactors = F, row.names = 1)
Summary_new = Summary[,c(1:3,8)]
write.csv(Summary_new,file="CD4_gene_subregion_DMS.csv",quote=F)
Count = data.frame(sum(Summary$Promoter),sum(Summary$Body),sum(Summary$UTR3),sum(Summary$Intergenic))
CD4 = read.csv("cg.comparison_aftercombatM.csv", header = T, stringsAsFactors = F)
CD4_DMS = CD4[which(CD4[,3] < 0.01),]
Gene = CD4_DMS$UCSC_RefGene_Group
Prop = c()
for (i in 1:length(Count)) {
  Prop = c(Prop,Count[i]/length(Gene))
}
Prop = t(as.data.frame(Prop))
Genefeature = data.frame(c("Promoter","Body","3'UTR","Intergenic"))
Geneprop = cbind(Prop,Genefeature)
Geneprop$Methylation_Sites = "Differential"
#All data
Summary_all = read.csv("Genefeature_aftercombatM.csv", header = T, stringsAsFactors = F)
Summary_all_new = Summary_all[,c(1:3,8)]
write.csv(Summary_all_new,file="CD4_gene_subregion_all.csv",quote=F)
Count_all = data.frame(sum(Summary_all$Promoter),sum(Summary_all$Body),sum(Summary_all$UTR3),sum(Summary_all$Intergenic))
anno_all = read.csv("anno_aftercombatM.csv", header = T, stringsAsFactors = F)
Gene_all = anno_all$UCSC_RefGene_Group
Prop_all = c()
for (i in 1:length(Count_all)) {
  Prop_all = c(Prop_all,Count_all[i]/length(Gene_all))
}
Prop_all = t(as.data.frame(Prop_all))
Genefeature = data.frame(c("Promoter","Body","3'UTR","Intergenic"))
Geneprop_all = cbind(Prop_all,Genefeature)
Geneprop_all$Methylation_Sites = "All"
colnames(Geneprop_all)= colnames(Geneprop)
Geneprop_merge = rbind(Geneprop,Geneprop_all)
colnames(Geneprop_merge)[1] = "Proportion"
colnames(Geneprop_merge)[2] = "Gene"
barplot = ggplot(Geneprop_merge, aes(x = Gene, y = Proportion, fill = Methylation_Sites)) + geom_bar(stat="identity", position="dodge") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 0.6, 0.1),limits = c(0,0.6)) +
  theme_classic()+
  theme(legend.position=c(0.86,0.90))+theme(legend.text = element_text(colour="black", size = 11))+theme(legend.background = element_rect(colour = "black"))+
  theme(axis.text.x=element_text(colour = "black", size=11))+ theme(axis.text.y=element_text(colour = "black", size=11))+
  theme(axis.title.x=element_blank()) +  theme(axis.title.y= element_text(colour = "black", size=12))+
  scale_color_npg()+scale_fill_npg()
ggsave(filename = "Genefeature_CD4_dva_aftercombatM_001_npg.tiff", plot = barplot, width = 7, height = 6, dpi = 500)
ggsave(filename = "Genefeature_CD4_dva_aftercombatM_001_npg.pdf", plot = barplot, width = 7, height = 6, dpi = 500)
