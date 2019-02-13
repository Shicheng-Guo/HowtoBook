cp ~/hpc/db/hg19/cpgSNP.hg19.bed ./
awk '{print "chr"$1,$2-1,$2,$3,$4,$5,$6,$7,$8}' OFS="\t" ~/hpc/db/hg19/cpgSNP.hg19.bed > cpgSNP.hg19.bed
awk '{print "chr"$1,$2-1,$2,$3}' OFS="\t" gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.bed > gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.hg19.bed
bedtools intersect -wo -a gnomad.genomes.r2.1.sites.rec.eQTL.set2.merge.vcf.hg19.bed -b cpgSNP.hg19.bed >  gnomad.genomes.eQTL.cpgSNP.hg19.bed
sort -u gnomad.genomes.eQTL.cpgSNP.hg19.bed > gnomad.genomes.eQTL.cpgSNP.uni.hg19.bed
wc -l gnomad.genomes.eQTL.cpgSNP.uni.hg19.bed
