time(bcftools view -T ../chr22.cpgSNP.bin.bed ../../chr22.1kg.phase3.v5a.vcf.gz -Oz -o chr22.cpgSNP.vcf.gz)
time(bcftools view -R ../chr22.cpgSNP.bin.bed ../../chr22.1kg.phase3.v5a.vcf.gz -Oz -o chr22.cpgSNP.R.vcf.gz)

method 1: 
time(bcftools view -T ../chr22.cpgSNP.bin.bed ../../chr22.1kg.phase3.v5a.vcf.gz -Oz -o chr22.cpgSNP.vcf.gz)
real    2m33.106s
user    2m32.641s
sys     0m0.293s

method 2:
bgzip ../chr22.cpgSNP.bin.bed
tabix -p bed ../chr22.cpgSNP.bin.bed.gz
time(bcftools view -R ../chr22.cpgSNP.bin.bed.gz ../../chr22.1kg.phase3.v5a.vcf.gz -Oz -o chr22.cpgSNP.vcf.gz)
real    9m4.696s
user    9m4.271s
sys     0m0.372s
