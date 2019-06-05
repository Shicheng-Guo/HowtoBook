```
for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr$i.vcf.bgz.tbi
done

for i in {1..22} X Y
do
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz.tbi
done


wget -O - "https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz" | gunzip -c | grep -E '(^#|\|(GSC|GSC2|FGF6|FSTL1)\|)' > genes.vcf
```
