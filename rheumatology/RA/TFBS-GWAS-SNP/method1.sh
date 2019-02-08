bedtools sort -i  GWAS-RA-792.R2.6.rsSNP.hg19.bed > GWAS-RA-792.R2.6.rsSNP.sort.hg19.bed
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.hg19.bed -b wgEncodeRegTfbsClusteredV3.hg19.bed | sort -u > GWAS-RA-792.R2.6.rsSNP.sort.tfbs.hg19.bed
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.hg19.bed -b wgEncodeRegDnaseClusteredV3.hg19.bed | sort -u > GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed 
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed  -b ~/hpc/db/hg19/CpGI_Shore.hg19.bed | sort -u >  GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.hg19.bed 
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed  -b ~/hpc/db/hg19/CpGI_Shelf.hg19.bed | sort -u >  GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.hg19.bed 
bedtools intersect -wa -a GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed  -b ~/hpc/db/hg19/CpGI.hg19.bed | sort -u >  GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.hg19.bed 

wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.hg19.bed
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.hg19.bed
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.hg19.bed 
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.hg19.bed 
wc -l GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.hg19.bed 
