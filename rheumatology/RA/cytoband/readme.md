Genetic Association SNPs from GWAS-Meta analysis (20190208)

You will find there are some non-related GWAS records come out if you use [grep](Grep_RA.GWAS_Catalog.md) to achieve RA-GWAS-SNPs. With manually check, we fix the number to [192 RA cytoband regions](RA-GWAS-Cytoband.hg19.bed) and then we collect all the missense/frame/stop functional SNPs within these 192 RA-cytoband regions, we found there are [19456 functional SNPs](gnomad.exomes.r2.1.sites.rec.RA-GWAS-Cytoband.hg19.vcf.bed). 

However, as Yukinori Okada mentioned in his [ARD 2018 paper](https://ard.bmj.com/content/early/2018/12/08/annrheumdis-2018-213678), Only[128 SNPs](../GWAS/GWAS-Meta-128-SNPs.20190208.vcf) have been identified in previous meta-analysis to GWAS-RA data. Steven told me this ~101 loci is not 101 SNPs, we should take them as 101 regions.
