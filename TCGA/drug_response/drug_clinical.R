library("TCGAbiolinks")
pid<-TCGAbiolinks:::getGDCprojects()$project_id
pid<-pid[grep("TCGA",pid)]
for(i in pid){
query <- GDCquery(project=i,data.category = "Clinical",file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query,"patient")
clinical.drug <- GDCprepare_clinic(query,"drug")
}
