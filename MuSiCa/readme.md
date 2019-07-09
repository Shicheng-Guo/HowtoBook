
```
if(!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if(!require(MutationalPatterns)) {BiocManager::install("MutationalPatterns")}
if(!require(VariantAnnotation)) {BiocManager::install("VariantAnnotation")}
if(!require(BSgenome.Hsapiens.UCSC.hg38)) {BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")}
if(!require(BSgenome.Hsapiens.UCSC.hg19)) {BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")}
if(!require(BSgenome.Hsapiens.1000genomes.hs37d5)) {BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")}

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(heatmaply)) install.packages("heatmaply")
if(!require(gplots)) install.packages("gplots")
if(!require(reshape2)) install.packages("reshape2")
if(!require(data.table)) install.packages("data.table")
if(!require(readxl)) install.packages("readxl")
if(!require(openxlsx)) install.packages("openxlsx")

if(!require(shiny)) install.packages("shiny")
if(!require(shinyBS)) install.packages("shinyBS")
if(!require(devtools)) install.packages("devtools")
if(!require(shinysky)) {library(devtools); devtools::install_github("AnalytixWare/ShinySky")}
if(!require(shinyjs)) install.packages("shinyjs")
if(!require(V8)) install.packages("V8")
if(!require(shinythemes)) install.packages("shinythemes")
if(!require(plotly)) install.packages("plotly")
if(!require(webshot)) install.packages("webshot")
webshot::install_phantomjs()
```
