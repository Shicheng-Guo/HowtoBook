
```
wget -O - "https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz" | gunzip -c | grep -E '(^#|\|(GSC|GSC2|FGF6|FSTL1)\|)' > genes.vcf
```
