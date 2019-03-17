
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/LungBrainMetastasis/bed")

cp LBM.MutationProfile.txt LBM.MutationProfile.test.txt
perl -p -i -e 's/3_prime_UTR_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/5_prime_UTR_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/coding_sequence_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/downstream_gene_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/frameshift_variant/frameshift_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/inframe_deletion/frameshift_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/intergenic_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/intron_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/mature_miRNA_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/missense_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/NMD_transcript_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/non_coding_transcript_exon_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/non_coding_transcript_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/protein_altering_variant/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/regulatory_region_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_acceptor_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_donor_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/splice_region_variant/splice_related_variants/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/start_lost/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/stop_gained/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/stop_lost/missense_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/synonymous_variant/Others/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/TF_binding_site_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt
perl -p -i -e 's/upstream_gene_variant/regulatory_region_variant/g' LBM.MutationProfile.test.txt

data<-read.table("LBM.MutationProfile.test.txt",head=T,sep="\t",row.names = 1)
data[1:3, 1:3]
target<-unique(c("TP53","KRAS","FAT4","STK11","EGFR","KMT2C","CHEK2P2","ERF","MIR4436A","BAGE2","BAGE3","BAGE4","BAGE5",
          "LOC102723769","MIR6077","FER1L4","FRG1BP","FRG1DP","LOC100996724","AHNAK2","FRG1BP","HMCN2","TPTE","BRSK1",
          "CR1L","CROCCP2","DUS3L","FRG1","MUC17","NOTCH2NL","SDHAP2","SSPO","ALPP","B4GALNT4","BAIAP3","BRSK2","MUC5B"))
input<-t(data[,c(na.omit(match(target,colnames(data))))])
input[1:3, 1:3]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  missense_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  regulatory_region_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  splice_related_variants = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#008000", col = NA))
  },
  frameshift_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "green", col = NA))
  },
  Others = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "purple", col = NA))
  }
  
)

col = c("missense_variant" = "red", "regulatory_region_variant" = "blue", "splice_related_variants" = "#008000","frameshift_variant"="green","Others"="purple")

oncoPrint(input, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col=col,
          remove_empty_columns = TRUE,
          heatmap_legend_param = list(title = "Alternations", at = c("missense_variant", "regulatory_region_variant", "splice_related_variants","frameshift_variant","Others"), 
                                      labels = c("missense_variant", "regulatory_region_variant", "splice_related_variants","frameshift_variant","Others")))


