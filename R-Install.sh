
source("http://bioconductor.org/biocLite.R")
biocLite("ggthemes")
biocLite("ggsci")
biocLite("dplyr")
biocLite("reshape")
biocLite("readr")
biocLite("readxl")

system("wget https://cran.fhcrc.org/src/contrib/ggthemes_3.2.0.tar.gz")
system("wget https://cran.r-project.org/src/contrib/ggsci_1.5.tar.gz")
system("wget https://cran.r-project.org/src/contrib/dplyr_0.5.0.tar.gz")
system("wget https://cran.r-project.org/src/contrib/reshape_0.8.6.tar.gz")
system("wget https://cran.r-project.org/src/contrib/hms_0.2.tar.gz")
system("wget https://cran.r-project.org/src/contrib/readr_1.0.0.tar.gz")
system("wget https://cran.r-project.org/src/contrib/readxl_0.1.1.tar.gz")

install.packages("ggsci_1.5.tar.gz")
install.packages("dplyr_0.5.0.tar.gz")
install.packages("reshape_0.8.6.tar.gz")
install.packages("readxl_0.1.1.tar.gz")
install.packages("hms_0.2.tar.gz")
install.packages("readr_1.0.0.tar.gz")


