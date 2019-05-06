
```
bedtools intersect -wao -b /home/guosa/hpc/project/pmrp/phase1/plink/RA/MCRI_RA_C2_Exom1_PvalueSort.hg19.bed -a /home/guosa/hpc/project/pmrp/phase2/RA/C2/MCRI_RA_C2_Exom2.hg19.bed > MCRI_RA_C2_Exom1_Exom2.hg19.bed


plink --meta-analysis /home/guosa/hpc/project/pmrp/phase1/plink/RA/PheTyp1_RA_C2.assoc /home/guosa/hpc/project/pmrp/phase2/RA/C2/PheTyp1_RA_C2.Exome1.assoc --out PheTyp1_RA_C2_Exom1_Exom2

/home/guosa/hpc/db/hg19/allSNP151.hg19.sort.bed


awk '{print "chr"$5,$6-1,$6,$1,$2,$3,$4,$5,$6,$7,$8}' OFS="\t" ExomeChip_SNPsInfo_hg19.txt | grep -v 'MapInfo'> ExomeChip_SNPsInfo_sort_hg19.txt


bedtools intersect -wao -a FinalRelease_QC_20140311_Team1_Marshfield.bim.bed -b /home/guosa/hpc/db/hg19/ExomeChip_SNPsInfo_sort_hg19.txt > FinalRelease_QC_20140311_Team1_Marshfield.bim.hg19.bed
bedtools intersect -wao -a FinalRelease_QC_20140311_Team1_Marshfield.bim.bed -b /home/guosa/hpc/db/hg19/allSNP151.hg19.sort.bed >  FinalRelease_QC_20140311_Team1_Marshfield.bim.hg19.bed

awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13}' OFS="\t" FinalRelease_QC_20140311_Team1_Marshfield.bim.hg19.bed | less -S 

plink --bfile FinalRelease_QC_20140311_Team1_Marshfield --allow-no-sex --pheno FinalRelease_QC_20140311_Team1_Marshfield.phen --pheno-name PheTyp1_RA_C2 --assoc counts --ci 0.95 --out ./RA/MCRI_RA_C2_Exom1
plink --bfile PMRP.PhaseII.Steven.Guo.RA.CEU.C2 --allow-no-sex --assoc counts --ci 0.95 --out MCRI_RA_C2_Exom2

plink --meta-analysis /home/guosa/hpc/project/pmrp/phase1/plink/RA/MCRI_RA_C2_Exom1.assoc /home/guosa/hpc/project/pmrp/phase2/RA/C2/MCRI_RA_C2_Exom2.assoc --out MCRI.PheTyp1_RA_C2_Exom1_Exom2


MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/project/pmrp/phase2/RA/C2")
library("Haplin")
png("qqplot.iron.png")
pQQ(na.omit(iron.SSc$P.Value), nlabs =nrow(na.omit(iron.SSc)), conf = 0.95)
dev.off()

awk '{print "chr"$1,$2-1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' OFS="\t" MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta.PvalueSort.hg19.bed | grep -v SNP > MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta.PvalueSort.hg19.uni.bed
bedtools intersect -wao -a MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta.PvalueSort.hg19.uni.bed -b ~/hpc/db/hg19/refGeneV2.hg19.bed | awk '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$19}' OFS="\t" > MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta.PvalueSort.hg19.RefGene.bed
head -n 200 MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta.PvalueSort.hg19.RefGene.bed | sort -u > MCRI.PheTyp1_RA_C2_Exom1_Exom2.meta.PvalueSort.hg19.RefGene.SigGH.txt
 
 ```
