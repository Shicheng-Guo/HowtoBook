How to download dbsnp151 vcf file in hg19

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.md5

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi

mv All_20180423.vcf.gz All_20180423.hg19.vcf.gz

mv All_20180423.vcf.gz.tbi All_20180423.hg19.vcf.gz.tbi

scp nu_guos@submit-1.chtc.wisc.edu:/home/nu_guos/All_2018* ./
