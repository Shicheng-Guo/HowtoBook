How to download vcf files (hg19 and hg38) for dbSNP152 

wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.bgz -O ~/hpc/db/hg19/dbSNP152.hg19.vcf.bgz
gunzip dbSNP152.hg19.vcf.bgz
wget https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.bgz -O ~/hpc/db/hg38/dbSNP152.hg38.vcf.bgz
gunzip dbSNP152.hg38.vcf.bgz
