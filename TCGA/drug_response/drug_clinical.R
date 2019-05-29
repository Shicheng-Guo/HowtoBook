library("TCGAbiolinks")
pid<-TCGAbiolinks:::getGDCprojects()$project_id
pid<-pid[grep("TCGA",pid)]

drug2csv<-function(clinical.drug){
  bcr_patient_barcode<-clinical.drug$bcr_patient_barcode
  therapy_types<-clinical.drug$therapy_types
  drug_name<-clinical.drug$drug_name
  measure_of_response<-clinical.drug$measure_of_response
  days_to_drug_therapy_start<-clinical.drug$days_to_drug_therapy_start
  days_to_drug_therapy_end<-clinical.drug$days_to_drug_therapy_end
  therapy_ongoing<-clinical.drug$therapy_ongoing
  new.clinical.drug<-data.frame(bcr_patient_barcode,therapy_types,drug_name,measure_of_response,days_to_drug_therapy_start,days_to_drug_therapy_end,therapy_ongoing)
  return(new.clinical.drug)
}

for(i in pid){
query <- GDCquery(project=i,data.category = "Clinical",file.type = "xml")
GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query,"drug")
drugResponse<-drug2csv(clinical.drug)
write.table(drugResponse,file=paste(i,".drugResponse.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
save(drugResponse,file=paste(i,".drugResponse.RData",sep=""))
}


