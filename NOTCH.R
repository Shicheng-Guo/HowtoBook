
ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi

cd /home/local/MFLDCLIN/guosa/hpc/db/1000Genome
plink --bfile chr6 --make-bed --keep /mnt/bigdata/Genetic/Projects/shg047/notch4/ChineseNonRel.txt --maf 0.05 --snps-only --chr 6 --from-bp 28601768 --to-bp 33472209  
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out notch4.input
plink --bfile notch4.input --r2 --matrix
20263
20363
awk '20263<=NR && NR <=20363' plink.ld > notch4.R2
cut -d' ' -f 20263-20363 plink.ld > Notch.R1.txt

32112 MB RAM detected; reserving 16056 MB for main workspace.
156119 out of 5024119 variants loaded from .bim file.
2504 people (0 males, 0 females, 2504 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink.nosex .
--keep: 144 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 144 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.99946.
156119 variants and 144 people pass filters and QC.
Note: No phenotypes present.
--make-bed to plink.bed + plink.bim + plink.fam ... done.
Pruned 134051 variants from chromosome 6, leaving 22068.
Pruning complete.  134051 of 156119 variants removed.
Marker lists written to plink.prune.in and plink.prune.out .

NOCH4 chr6:32162620-32191844

A set of oligonucleotide probes was used to capture the ∼5-Mbp MHC region from sheared genomic DNA libraries. The targeted region is located at Chromosome 6p21 and encompasses all genes from GPX5 to ZBTB9 (Fig. 1A)

GPX5 chr6:28,601,768-28,609,923
ZBTB9 chr6:33,469,263-33,472,209

28601768
33472209


In a simulation study,
Kruglyak (1999) found that LD was unlikely to extend
beyond an average of 3 kb in general populations and
in most isolated populations, so that >=500,000 SNPs
would be required for whole-genome association studies.
On the other hand, Reich et al. (2001) showed that
LD in a U.S. population of northern European descent
could extend 60 kb for common alleles, so that only
50,000 SNPs would be needed in these populations.
