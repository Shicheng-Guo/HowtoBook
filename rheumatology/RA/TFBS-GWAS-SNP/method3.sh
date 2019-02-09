awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' GWAS-Meta-128-SNPs.20190208.vcf | sort > GWAS-Meta-128-SNPs.20190208.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/rheumatology/RA/TFBS-GWAS-SNP/RA-GWAS-Functional-SNPs-Final.186.Snp.20190208.bed
bedtools intersect -a RA-GWAS-Functional-SNPs-Final.186.Snp.20190208.bed -b GWAS-Meta-128-SNPs.20190208.bed > RA-GWAS-Functional-SNPs-Final.186.Snp.GWAS-Meta-129.overlap.20190208.bed
