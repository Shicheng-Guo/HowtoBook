

miRNA target prediction by [targetScan](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi) can be download with [default mode](http://www.targetscan.org/vert_72/vert_72_data_download/Predicted_Target_Locations.default_predictions.hg19.bed.zip) and [full mode](http://www.targetscan.org/vert_72/vert_72_data_download/Conserved_Site_Context_Scores.txt.zip).



miRNA-SNP: http://bioinfo.life.hust.edu.cn/miRNASNP2/download.php

* miRNA_gain_by_SNPs_in_gene_3utr.txt
* miRNA_loss_by_SNPs_in_gene_3utr.txt
* miRNA_targets_loss_by_SNPs_in_seed_regions.txt
* miRNA_targets_gain_by_SNPs_in_seed_regions.txt
```
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_gain_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_targets_loss_by_SNPs_in_seed_regions.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_gain_by_SNPs_in_gene_3utr.txt
wget http://bioinfo.life.hust.edu.cn/miRNASNP2/download/miRNA_loss_by_SNPs_in_gene_3utr.txt
```
Prepare shared SNPs
```
setwd("/gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/All_Target_Locations.hg19.bed")
d1<-read.table("/gpfs/home/guosa/hpc/rheumatology/RA/KEGG_90_RA_GeneList.txt")
d2<-read.table("/gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/Design/Target.miRNA.mature.bed")
data<-c()
file=list.files(pattern="*.txt")
for(i in 1:length(file)){
print(i)
d3<-read.table(file[i])
d4<-d3[d3[,4]%in%d1[,1],]
d5<-d4[d4[,5]%in%d2[,8],]
data<-rbind(data,d5)
}
```
