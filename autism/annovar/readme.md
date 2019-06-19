Step 1. download dbnsfp30a -> dbnsfp40a (hg19, set ppn=8)
```
annotate_variation.pl -downdb -webfrom annovar -build hg19 dbnsfp30a   ~/hpc/tools/annovar/humandb/
annotate_variation.pl -downdb -webfrom annovar -build hg19 dbnsfp40a   ~/hpc/tools/annovar/humandb/
annotate_variation.pl --buildver hg19 --downdb seq /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_seq
retrieve_seq_from_fasta.pl /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_ensGene.txt -seqdir /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_seq -format ensGene -outfile /gpfs/home/guosa/hpc/tools/annovar/humandb/hg19_ensGeneMrna.fa
```
Step 2. annotation SNPs with annovar (dbnsfp35c,hg19, set ppn=8)
```
mkdir annovar
mkdir temp
for i in {1..22} 
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=8 >> chr$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> chr$i.job
echo \#PBS -m abe  >> chr$i.job
echo \#PBS -o $(pwd)/temp/ >>chr$i.job
echo \#PBS -e $(pwd)/temp/ >>chr$i.job
echo cd $(pwd) >> chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq All_samples_Exome_QC.chr$i.vcf.vcf.gz  \> ./annovar/All_samples_Exome_QC.chr$i.avinput >> chr$i.job
echo table_annovar.pl ./annovar/All_samples_Exome_QC.chr$i.avinput /gpfs/home/guosa/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 --csvout -out ./annovar/chr$i -remove -protocol refGene,dbnsfp33a,gwasCatalog,wgRna,targetScanS,tfbsConsSites -operation gx,f,r,r,r,r -nastring . -otherinfo -polish -xref /gpfs/home/guosa/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done
```
Step 3. filter delterious SNPs with dbnsfp35c (6 out of 11 as deleterious)
```
RiskAAC<-function(data){
  risk1<-grep("D",data$SIFT_pred)
  risk2<-grep("D|P",data$Polyphen2_HDIV_pred)
  risk3<-grep("D|P",data$Polyphen2_HVAR_pred)
  risk4<-grep("D",data$LRT_pred)
  risk5<-grep("D",data$MutationTaster_pred)
  risk6<-grep("H|M",data$MutationAssessor_pred)
  risk7<-grep("D",data$FATHMM_pred)
  risk8<-grep("D",data$PROVEAN_pred)
  risk9<-grep("D",data$MetaSVM_pred)
  risk10<-grep("D",data$MetaLR_pred)
  risk11<-grep("D",data$fathmm.MKL_coding_pred)
  risk12<-grep("D",data$M.CAP_pred)
  rlt=c(risk1,risk2,risk3,risk4,risk5,risk6,risk7,risk8,risk9,risk10,risk11,risk12)
  return(rlt)
}

file=list.files(pattern="*.csv")
rlt<-c()
for(i in 1:length(file)){
  data<-read.csv(file[i])
  Riskloci<-names(which(table(RiskAAC(data))>=6))
  rlt<-rbind(rlt,data[Riskloci,])
  print(i)
}
dim(rlt)
write.table(rlt,file=paste("autism","anno.10748.txt",sep="."),col.names=F,row.names=F,quote=F)
```
