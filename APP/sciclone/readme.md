Install sciClone in R-3.6.0 in Windows
```
BiocManager::install("devtools")
BiocManager::install("IRanges")
BiocManager::install("digest")
BiocManager::install("IRanges")
install.packages("https://cran.rstudio.com/bin/windows/contrib/3.6/processx_3.3.1.zip", repos = NULL)
install.packages("https://cran.rstudio.com/bin/windows/contrib/3.6/callr_3.2.0.zip", repos = NULL)

library('processx')
library("callr")
library("devtools")
library("IRanges")
install_github("genome/bmm")
install_github("genome/sciClone")
```
