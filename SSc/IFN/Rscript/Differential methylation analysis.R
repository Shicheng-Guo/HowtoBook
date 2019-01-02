#Preprocessing of idat files - example
  working_dir = ""
  sample_annotation = ""
  library(RnBeads)
  library(doParallel)
  #parameter setting 
  logger.start(fname=NA)
  num.cores=8;
  parallel.setup(num.cores);
  print.gtable=grid.draw
  data.dir = working_dir
  idat.dir = file.path(data.dir, "idat_CD8")
  sample.annotation = sample_annotation
  analysis.dir = paste(data.dir,"analysis",sep="/")
  report.dir = file.path(analysis.dir, "reports")
  # data import 
  data.source = c(idat.dir, sample.annotation);
  result = rnb.run.import(data.source = data.source,data.type="infinium.idat.dir", dir.reports = report.dir)
  print("Finished data importing")
  # data preprocessing including filtering and normalization
  rnb.options(filtering.sex.chromosomes.removal=TRUE,normalization.method = "bmiq",filtering.greedycut.pvalue.threshold=0.05,filtering.missing.value.quantile=0,
              filtering.deviation.threshold=0,filtering.snp="any",filtering.cross.reactive = TRUE);
  rnb.set.unfiltered = result$rnb.set
  result = rnb.run.preprocessing(rnb.set.unfiltered,dir.reports = report.dir);
  rnb.set = result$rnb.set;
  methydata=meth(rnb.set,row.names=TRUE);
  save(methydata,file="GD_CD8_data_beforecombat.RData")
  
#Batch normalization for integrated methylation data - example
 load("CD8_GS_data_beforecombat.RData")
 library(lumi)
 CD8_GS_data_beforecombat_M = beta2m(CD8_GS_data_beforecombat)
 matrix_c = data.matrix(CD8_GS_data_beforecombat_M)
 clinical = read.csv("CD8_GS_clinical.csv", header = T)
 library(sva)
 CD8_GS_data_aftercombat_M = ComBat(matrix_c,as.factor(clinical$Disease_Type))
 CD8_GS_data_aftercombatM = m2beta(CD8_GS_data_aftercombat_M)
 write.csv(CD8_GS_data_aftercombatM, file = "CD8_GS_data_aftercombatM.csv", quote = F)
 save(CD8_GS_data_aftercombatM,file="CD8_GS_data_aftercombatM.RData")

#Evaluate methylation difference between case and control - example
diffcomp = function(working_directory="",data_matrix="CD4_RGSS_data_aftercombatM.RData",
                    clinical_file="CD4_RGSS_clinical.csv", case="Case", control="Control"){
  setwd(working_directory);
  load(data_matrix);
  clinical = read.csv(clinical_file, header = T);
  seq_case = ifelse(clinical$Sample_Group==case,TRUE,FALSE);
  seq_control = ifelse(clinical$Sample_Group==control,TRUE,FALSE);
  meandiff = apply(CD4_RGSS_data_aftercombatM,1,function(x){return(abs(mean(x[seq_case],na.rm = T) - mean(x[seq_control],na.rm = T) ))})
  foldchange =apply(CD4_RGSS_data_aftercombatM,1,function(x){return( mean(x[seq_case],na.rm = T) / mean(x[seq_control],na.rm = T) )});
  
  #creating variates: Disease_Type, Age and Gender
  x = ifelse(clinical$Sample_Group == case ,1,0);
  Age = clinical$Age;
  Gender = clinical$Gender;
  Eth = clinical$Eth;
  
  pvalue =apply(CD4_RGSS_data_aftercombatM, 1, function(y){return( summary(lm(y~as.factor(x)+Age+as.factor(Gender)+as.factor(Eth)))$coeff[2,4])})
  mean_case =apply(CD4_RGSS_data_aftercombatM, 1, function(x) {return( mean(x[seq_case],na.rm = T ))});
  mean_control = apply(CD4_RGSS_data_aftercombatM, 1, function(x) {return( mean(x[seq_control],na.rm = T ))});
  pvalue.fdr = p.adjust(pvalue,method = "fdr");
  cg.comparison = data.frame(pvalue,pvalue.fdr,mean_case,mean_control,meandiff,foldchange);
  cg.comparison = cg.comparison[order(cg.comparison[,1]),];
  library(readr)
  anno = read_tsv("HM450KAnnotation.txt")
  seq = match(rownames(cg.comparison),anno$ID);
  anno_top = anno[seq,c(11,12,21,23,25)];
  cg.comparison_anno = data.frame(cg.comparison, anno_top);
  write.csv(cg.comparison_anno,file="cg.comparison_aftercombatM.csv",quote=F);
  
}