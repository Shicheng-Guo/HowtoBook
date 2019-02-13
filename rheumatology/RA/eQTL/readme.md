

In the first step, we collected all the eQTL from whole blood, liver, colon, stomach [etql set1 code](eqtl.set1.sh). 
* [21098](eQTL.hg19.bed) eQTL SNPs were collected from the GTEx project with FDR<0.05
* [14285](gnomad.genomes.r2.1.sites.rec.eQTL.merge.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian
* [772](gnomad.exomes.r2.1.sites.rec.eQTL.hg19.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian

However, when I check it again, I found I lost the lung eqtl data. After I add lung eqtl [eqtl set2 code](eqtl.set2.sh)
* [30517](eQTL.set2.hg19.bed) eQTL SNPs were collected from the GTEx project with FDR<0.05
* [20667](gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian
* [772](gnomad.exomes.r2.1.sites.rec.eQTL.set2.hg19.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian

I made a brief checking to these eqtl candidates.

* in 30517 is CpG-SNP
* in 30517 is CpG-SNP
* in 30517 is CpG-SNP

