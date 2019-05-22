if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")

library("TCGAbiolinks")
library("data.table")

followup2csv<-function(clinical.fellowup){
  bcr_patient_barcode=clinical.fellowup$bcr_patient_barcode
  lost_follow_up=clinical.fellowup$lost_follow_up
  vital_status=clinical.fellowup$vital_status
  primary_therapy_outcome_success=clinical.fellowup$primary_therapy_outcome_success
  days_to_death=clinical.fellowup$days_to_death
  new_tumor_event_after_initial_treatment=clinical.fellowup$new_tumor_event_after_initial_treatment
  days_to_new_tumor_event_after_initial_treatment=clinical.fellowup$days_to_new_tumor_event_after_initial_treatment
  person_neoplasm_cancer_status=clinical.fellowup$person_neoplasm_cancer_status
  karnofsky_performance_score=clinical.fellowup$karnofsky_performance_score
  new.clinical.fellowup<-data.frame(bcr_patient_barcode,lost_follow_up,
                                    vital_status,primary_therapy_outcome_success,days_to_death,
                                    new_tumor_event_after_initial_treatment, 
                                    days_to_new_tumor_event_after_initial_treatment,person_neoplasm_cancer_status,karnofsky_performance_score)
  return(new.clinical.fellowup)
}
patient2csv<-function(clinical.patient){
  bcr_patient_barcode<-clinical.patient$bcr_patient_barcode
  gender<-clinical.patient$gender
  tobacco_smoking_history<-clinical.patient$tobacco_smoking_history
  number_pack_years_smoked<-clinical.patient$number_pack_years_smoked
  year_of_tobacco_smoking_onset<-clinical.patient$year_of_tobacco_smoking_onset
  days_to_last_followup<-clinical.patient$days_to_last_followup
  days_to_death<-clinical.patient$days_to_death
  person_neoplasm_cancer_status<-clinical.patient$person_neoplasm_cancer_status
  age_at_initial_pathologic_diagnosis<-clinical.patient$age_at_initial_pathologic_diagnosis
  anatomic_neoplasm_subdivision<-clinical.patient$anatomic_neoplasm_subdivision
  diagnosis<-clinical.patient$diagnosis
  stage_event_tnm_categories<-clinical.patient$stage_event_tnm_categories
  stage_event_pathologic_stage<-clinical.patient$stage_event_pathologic_stage
  new.clinical.patient<-data.frame(bcr_patient_barcode,gender,age_at_initial_pathologic_diagnosis,tobacco_smoking_history,number_pack_years_smoked,
                                   year_of_tobacco_smoking_onset,
                                   days_to_last_followup, days_to_death,person_neoplasm_cancer_status,
                                   anatomic_neoplasm_subdivision,diagnosis,stage_event_tnm_categories,stage_event_pathologic_stage)
  return(new.clinical.patient)
}
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
radiation2csv<-function(clinical.radiation){
  bcr_patient_barcode<-clinical.radiation$bcr_patient_barcode
  radiation_type<-clinical.radiation$radiation_type
  radiation_dosage<-clinical.radiation$radiation_dosage
  radiation_treatment_ongoing<-clinical.radiation$radiation_treatment_ongoing
  measure_of_response<-clinical.radiation$measure_of_response
  days_to_radiation_therapy_start<-clinical.radiation$days_to_radiation_therapy_start
  days_to_radiation_therapy_end<-clinical.radiation$days_to_radiation_therapy_end
  measure_of_response<-clinical.radiation$measure_of_response
  new.clinical.radiation<-data.frame(bcr_patient_barcode,radiation_type,radiation_dosage,radiation_treatment_ongoing,days_to_radiation_therapy_start,days_to_radiation_therapy_end,measure_of_response)
  return(new.clinical.radiation)
}
project="TCGA-OV"
# download all xml files into current directory
query <- GDCquery(project =project,data.category = "Clinical", file.type = "xml")
GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
drug<-unique(drug2csv(clinical.drug))
head(drug)
save(drug,file=paste(project,"drug.RData",sep="."))
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
radiation<-unique(radiation2csv(clinical.radiation))
head(radiation)
save(radiation,file=paste(project,"radiation.RData",sep="."))
clinical.fellowup <- GDCprepare_clinic(query, clinical.info = "follow_up")
followup<-unique(followup2csv(clinical.fellowup))
head(followup)
save(followup,file=paste(project,"followup.RData",sep="."))
clinical.patient <- GDCprepare_clinic(query, clinical.info = "patient")
patient<-unique(patient2csv(clinical.patient))
head(patient)
save(patient,file=paste(project,"patient.RData",sep="."))
clinical.newtumorevent <- GDCprepare_clinic(query, clinical.info = "new_tumor_event")
newtumorevent<-unique(clinical.newtumorevent)
head(newtumorevent)
save(newtumorevent,file=paste(project,"newtumorevent.RData",sep="."))
clinical.stageevent <- GDCprepare_clinic(query, clinical.info = "stage_event")
stageevent<-unique(clinical.stageevent)
head(stageevent)
save(stageevent,file=paste(project,"stageevent.RData",sep="."))

