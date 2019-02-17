
Print Asian allele frequency
```
bcftools query  -f '%CHROM\t%POS\t%ID\t%INFO/controls_AF_eas\n' gnomad.exomes.r2.1.sites.chr1.rec.refGene.sort.rmdup.biallelic.vcf.bgz | awk '{print "chr"$1,$2-1,$2,$3,$4}' OFS="\t"
```
