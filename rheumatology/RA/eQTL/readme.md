

In the first step, we collected all the eQTL from whole blood, liver, colon, stomach [etql set1 code](eqtl.set1.sh). 
* [21098](eQTL.hg19.bed) eQTL SNPs were collected from the GTEx project with FDR<0.05
* [14285](gnomad.genomes.r2.1.sites.rec.eQTL.merge.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian
* [772](gnomad.exomes.r2.1.sites.rec.eQTL.hg19.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian

However, when I check it again, I found I lost the lung eqtl data. After I add lung eqtl [eqtl set2 code](eqtl.set2.sh)
* [30517](eQTL.set2.hg19.bed) eQTL SNPs were collected from the GTEx project with FDR<0.05
* [20667](gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian
* [1062](gnomad.exomes.r2.1.sites.rec.eQTL.set2.hg19.vcf.bed) SNPs identified by Gnomad.Genomes data with MAF>0.1% in East Asian

I made a brief checking to these eqtl candidates.[eqtl set3 code](eqtl.set3.sh)

* [6876](gnomad.genomes.eQTL.cpgSNP.uni.hg19.bed) in 20667 is belong to CpG-SNP  (Ratio=33%)
* [1199] SNPs were "Conflicting CpG-SNP" such as rs999941, rs9974367, they have pattern [C[C/G]G](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?do_not_redirect&rs=rs9974367)
* [5570](gnomad.genomes.eQTL.cpgSNP.uni.non-flict-CpG-SNP.hg19.bed) in 20667 SNPs are non-conflicted CpG-SNP eQTL and were added to ASA array. 


