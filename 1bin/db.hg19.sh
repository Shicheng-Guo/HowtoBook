How to start your new computational lab quickly

# 1. download dbsnp151,dbsnp152 vcf
# 2. 

# download dbsnp152 for hg19 and hg38
wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.bgz -O ~/hpc/db/hg19/dbSNP152.hg19.vcf.bgz
tabix -p vcf dbSNP152.hg19.vcf.bgz
zcat dbSNP152.hg19.vcf.bgz > dbSNP152.hg19.vcf
wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.bgz -O ~/hpc/db/hg38/dbSNP152.hg38.vcf.bgz
tabix -p vcf dbSNP152.hg38.vcf.bgz
zcat dbSNP152.hg19.vcf.bgz > dbSNP152.hg38.vcf.bgz

