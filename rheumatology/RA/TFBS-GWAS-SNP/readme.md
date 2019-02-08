DNA Methylation in GWAS Significant Rheumatoid Arthritis Associated Regions. 

We download 791 GWAS-Significant SNPs from GWAS Catalog and we collected all the linked SNPs with R2>0.6 in Asian Population. Totally, we identified 21,079 SNPs with above method. 

[Method 1](method1.sh):

* RA-LD-SNP (21079) -> TFBS (3766) -> DNase (2054) -> CpGisland (129) -> [129 SNPs](GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.hg19.bed) -> eQTL([19 SNP](GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI.eQTL.hg19.bed))
* RA-LD-SNP (21079) -> TFBS (3766) -> DNase (2054) -> CpG-Shore (316)-> [316 SNPs](GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.hg19.bed) -> eQTL([43 SNP](GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shore.eQTL.hg19.bed))
* RA-LD-SNP (21079) -> TFBS (3766) -> DNase (2054) -> CpG-Shelf (202)-> [202 SNPs](GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.hg19.bed) -> eQTL([39 SNP](GWAS-RA-792.R2.6.rsSNP.sort.tfbs.Dnase.CpGI_Shelf.eQTL.hg19.bed))

Summary (1): 10 multi-hits SNPs and [82 unique SNP](RA-GWAS-Functional-SNPs-Final.82.Snp.20190208.bed) and [55 genomic regions](RA-GWAS-Functional-SNPs-Final.82.Snp.20190208.sort.merge.bed) were occured in the final SNPs list. 

[Method 2](method2.md): 

* S1: R2>0.6 -> eQTL(PRA) -> TFBS -> DNase -> CpG-Island -> [7 Genomic Regions](S1-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed) and [7 SNPs](S1-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.bed)
* S2: R2>0.6 -> eQTL(Full) -> TFBS -> DNase -> CpG-Island -> [17 Genomic Regions](S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed) and [19 SNPs](S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.sort.bed)
* S3: R2>0.6 -> exom + missense + stop + frame + MAF_Asian>0.1% -> [104 SNPs](gnomad.exomes.r2.1.sites.rec.GWAS-RA-792.R2.6.rsSNP.input.hg19.vcf.bed)

Summary (1 & 2): 36 multi-hits SNPs and [186 unique SNP](RA-GWAS-Functional-SNPs-Final.186.Snp.20190208.bed) and [124 genomic regions](RA-GWAS-Functional-SNPs-Final.186.Snp.20190208.sort.merge.hg19.bed) were occured in the final SNPs list. 

[Method 3](method3.md): 
* S1: R2>0.6 -> eQTL(PRA, [206 SNPs](GWAS-RA-792.R2.6.rsSNP.PRA.eQTL.hg19.bed)) -> TFBS -> DNase -> CpG-Island -> [7 Genomic Regions](S1-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed) and [7 SNPs](S1-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.bed)
* S2: R2>0.6 -> eQTL(Full,[SNPs](GWAS-RA-792.R2.6.rsSNP.FullRA.eQTL.hg19.bed)) -> TFBS -> DNase -> CpG-Island -> [17 Genomic Regions](S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.merge.sort.bed) and [19 SNPs](S2-GWAS-RA-R2.6.eQTL.tfbs.DNase.CpGI.hg19.sort.bed)

