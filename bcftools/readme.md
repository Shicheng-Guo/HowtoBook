#### bcftools — utilities for variant calling and manipulating VCFs and BCFs.
#### version: 1.9- 2019-01-31

How to select bi-allelic SNPs from VCF files. (-d and -m cannot be used at the same time)
sort vcf, recover multiallelic, rm duplication, view -m2 -M2 -v snps
```
bcftools sort check.temp.17.vcf > check.17.vcf
bcftools norm  -m + check.17.vcf > check.17.recover.vcf
bcftools norm -d both check.17.recover.vcf > check.17.recover.2.vcf
bcftools view -m2 -M2 -v snps check.17.recover.2.vcf > check.17.trim.vcf
grep rs1799966 check.17.trim.vcf | less -S 
```
Time difference between different output format: u > v > b > z
```
time(bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.sort.bgz -Ob -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz)
Lines   total/split/realigned/skipped:  416866/0/0/0

real    2m23.143s
user    2m19.562s
sys     0m2.836s

time(bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.sort.bgz -Ou -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz)
Lines   total/split/realigned/skipped:  416866/0/0/0

real    1m16.727s
user    0m57.623s
sys     0m2.957s

time(bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.sort.bgz -Oz -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz)
Lines   total/split/realigned/skipped:  416866/0/0/0

real    5m32.718s
user    5m27.886s
sys     0m4.081s

time(bcftools norm -m + gnomad.exomes.r2.1.sites.chr22.vcf.sort.bgz -Ov -o gnomad.exomes.r2.1.sites.chr22.rec.vcf.bgz)
Lines   total/split/realigned/skipped:  416866/0/0/0

real    1m56.530s
user    1m26.564s
sys     0m6.377s
```
