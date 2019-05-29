query <- GDCquery(data.category = "Clinical",file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query,"patient")
clinical.drug <- GDCprepare_clinic(query,"drug")
