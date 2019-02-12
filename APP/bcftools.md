filter SNPs with rsID with bcftools
* bcftools view -i 'ID=@InnateDB.UTR3.rsid.txt' ~/hpc/db/hg19/All_20180423.vcf | grep -v '#' | awk '{print $1,$2,$3,$4,$5}' OFS="\t"
