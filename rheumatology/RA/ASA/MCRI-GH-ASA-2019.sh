
awk '{print "chr"$2,$3-1,$3,$4}' OFS="\t"  ASACHIA_rsID.txt | grep -v 'CHROM' | sort -u | bedtools sort > ASA.hg19.bed

bedtools intersect -v -a 15446.MRCI.ASA.eQTL.hg19.MAF0.001.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.genomes.r2.1.sites.rec.GWASCatalog.ASA.merge.vcf.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 4678.UTR3miRNAsNP.EAS.ImmuGene.MAF0.01.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 1169.UTR3miRNAsNP.EAS.MAF0.01.WGSC.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a GWAS-Meta-128-SNPs.20190208.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 48-SNPs-Genetic-Risk-Score-EUR.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.HLA.hg19.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 1375.gnomad.exomes.r2.1.sites.rec.RA-GWAS-Cytoband.hg19.vcf.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a GWAS-immnue-3325_SNP.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.40KGWASAID.merge.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.InnateDB.merge.hg19.bed -b ASA.hg19.bed | sort -u| awk '{print $1,$2,$3,$4}' OFS="\t" > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a cpgSNPisland.AID.GWAS.SNP.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a gnomad.exomes.r2.1.sites.rec.RAcandidate.hg19.vcf.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 63.miRNA.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 47_NG2012.SupplementaryTable5.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 145.UTR3miRNAsNP.EAS.MAF0.01.VIP.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a HLA-TagSnp.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a 1345.UTR3miRNAsNP.EAS.MAF0.01.EpiGeneC.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed

bedtools intersect -v -a miRNASeedSNP.hg19.bed -b ASA.hg19.bed | awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > ASA-GH-Guo.A1.bed
wc -l ASA-GH-Guo.A1.bed
cat ASA-GH-Guo.A1.bed >> ASA.hg19.bed


grep -v 'CHROM'  ASA.hg19.bed | sort -u > ASA.hg19.P2.bed
wc -l ASA.hg19.P2.bed
awk '{print "chr"$2,$3-1,$3,$4}' OFS="\t" ASACHIA_rsID.txt | grep -v CHROM | sort -u > ASA.hg19.P1.bed
wc -l ASA.hg19.P1.bed
bedtools intersect -v -a ASA.hg19.P2.bed  -b ASA.hg19.P1.bed| awk '{print $1,$2,$3,$4}' OFS="\t" | sort -u > MRCI.hg19.bed
wc -l MRCI.hg19.bed

