* DNA methylation (mh450) for UCEC in TCGA
```
setwd("/home/hgc/mxiong/ucec/MH450")
files=list.files(pattern="*8.txt.dat$",recursive = T)
methdata<-c()
for(i in 1:length(files)){
  temp<-read.table(files[i],head=T,sep="\t",row.names = 1)
  methdata<-cbind(methdata,temp[,1])
  print(i)
}
colnames(methdata)<-files
rownames(methdata)<-rownames(temp)
methdata[1:5,1:5]
save(methdata,file="TCGA.UCEC.mh450.RData")


setwd("/home/hgc/mxiong/ucec/miRNA")
files=list.files(pattern="*/*mirnas.quantification.txt$",recursive = T)
methdata<-c()
for(i in 1:length(files)){
  temp<-read.table(files[i],head=T,sep="\t",row.names = 1)
  methdata<-cbind(methdata,temp[,2])
  print(i)
}
colnames(methdata)<-files
rownames(methdata)<-rownames(temp)
methdata[1:5,1:5]
miRNA<-methdata
save(miRNA,file="TCGA.UCEC.miRNA.RData")

manifest="gdc_manifest.tcga.ucec.miRNA_BCGSC_miRNA_Profiling.2019-05-23.txt"
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste0("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system("curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > barcode.txt")

miRinfo<-read.table("miRNA_barcode.txt",head=T,sep="\t")
filename<-unlist(lapply(colnames(miRNA),function(x) unlist(strsplit(x,"/"))[2]))
newfilename<-unlist(sam[match(filename,sam$file_name),5])
colnames(miRNA)<-newfilename
save(miRNA,file="")
```
