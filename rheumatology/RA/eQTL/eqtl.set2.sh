/gpfs/home/guosa/hpc/rheumatology/RA/ASA/eqtl/SNP
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
tar xzvf GTEx_Analysis_v7_eQTL.tar.gz
cd GTEx_Analysis_v7_eQTL

qval_threshold=0.05
data1<-subset(read.table("Whole_Blood.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data2<-subset(read.table("Liver.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data3<-subset(read.table("Small_Intestine_Terminal_Ileum.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data4<-subset(read.table("Stomach.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
data5<-subset(read.table("Lung.v7.egenes.txt",head=T,sep="\t"),qval<qval_threshold)
eqtl<-c(as.character(data1[,19]),as.character(data2[,19]),as.character(data3[,19]),as.character(data4[,19]),as.character(data5[,19]))
length(table(eqtl))
eqtl.snp<-names(table(eqtl))
write.table(eqtl.snp,file="eqtl.snp.txt",sep="\t",quote=F,col.names=F,row.names=F)
