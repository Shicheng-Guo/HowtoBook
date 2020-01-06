dir.make("qcRNAseq")

# only if you don't have these packages already:
# install.packages("stringr")
# install.packages("gap")
# install.packages("sqldf")
# install.packages("reshape")
library(stringr) # makes string splitting, indexing, regexes easier - R's base functions are quite bad at this
library(gap) # best library I've found for qq plots
library(sqldf) # for me, JOIN & GROUP BY are easier to use than R's match/merge and aggregate
library(reshape) # for converting between matrices and relational models
options(stringsAsFactors=FALSE) # so that strings stay as strings in all read.table calls
