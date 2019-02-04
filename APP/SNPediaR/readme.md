```
BiocManager::install("SNPediaR", version = "3.8")
install.packages("devtools")
library(devtools)
install_github("genometra/SNPediaR")
library (SNPediaR)

res <- getCategoryElements(category = "Is_a_medical_condition")
head (res)
grep ('cancer', res, value = TRUE)

res <- getCategoryElements(category = "Is_a_medicine")
head (res)
grep ('cancer', res, value = TRUE)

res <- getCategoryElements(category = "Is_a_topic")
head (res)
grep ('cancer', res, value = TRUE)

res <- getCategoryElements(category = "Is_a_snp")
head (res)

res <- getCategoryElements(category = "In_dbSNP")
head (res)

res <- getCategoryElements(category = "Is_a_genotype")
head (res)
```
