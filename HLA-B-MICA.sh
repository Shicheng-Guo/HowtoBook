hg19 genomic region for HLA-B and MICA
chr6:31,232,075-31,391,038

cd /home/local/MFLDCLIN/guosa/hpc/db/1000Genome
plink --bfile chr6 --make-bed --keep /mnt/bigdata/Genetic/Projects/shg047/notch4/ChineseNonRel.txt --maf 0.05 --snps-only --chr 6 --from-bp 28601768 --to-bp 33472209  
awk '{print $2}' plink.bim > mysnps.txt
plink --bfile plink --list-all --show-tags mysnps.txt
plink --bfile plink --extract plink.tags --make-bed --recode --tab --out notch4.input
plink --bfile notch4.input --r2 --matrix
