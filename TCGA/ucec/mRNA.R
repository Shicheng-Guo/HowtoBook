files=list.files(pattern="*FPKM-UQ.txt$",recursive = T)
methdata<-c()
for(i in 1:length(files)){
  temp<-read.table(files[i],head=T,sep="\t",row.names = 1)
  methdata<-cbind(methdata,temp[,1])
  print(i)
}
colnames(methdata)<-files
rownames(methdata)<-rownames(temp)
methdata[1:5,1:5]
RNAseq<-methdata
save(RNAseq,file="TCGA.UCEC.RNAseq_FPKM.RData")

manifest="gdc_manifest.tcga.ucec.HTSeq-FPKM-UQ.2019-05-23.txt"
x=read.table(manifest,header = T)
manifest_length= nrow(x)
id= toString(sprintf('"%s"', x$id))
Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste0("\"",manifest_length, "\"", "}")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload.txt",quote=F,col.names=F,row.names=F)
system("curl --request POST --header \"Content-Type: application/json\" --data @Payload.txt \"https://api.gdc.cancer.gov/files\" > barcode.txt")

miRinfo<-read.table("FPKM_barcode.txt",head=T,sep="\t")
newfilename<-unlist(miRinfo[match(colnames(RNAseq),miRinfo$file_name),5])
colnames(RNAseq)<-newfilename
save(RNAseq,file="../TCGA.UCEC.RNAseq_FPKM.RData")
