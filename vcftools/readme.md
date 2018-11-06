Error: rs4565649 is monomorphic across all valid observations.

How to merge two vcf files: it is still the best way to merge vcf with plink
```
plink --bfile chr22.chip --list-duplicate-vars 
awk '{print $4}' plink.dupvar | grep -v ID > plink.dupvar.id 
plink --bfile chr22.chip --exclude plink.dupvar.id --make-bed --out chr22.chip.rmdup
plink --bfile chr22.imputation --list-duplicate-vars 
awk '{print $4}' plink.dupvar | grep -v ID > plink.dupvar.id 
plink --bfile chr22.imputation --exclude plink.dupvar.id --make-bed --out chr22.imputation.rmdup
plink --bfile chr22.imputation.rmdup --bmerge chr22.chip.rmdup --make-bed --out merge
plink --bfile chr22.chip.rmdup --flip merge-merge.missnp --make-bed --out chr22.chip.rmdup.flip
plink --bfile chr22.imputation.rmdup --bmerge chr22.chip.rmdup.flip --make-bed --out merge
plink --bfile chr22.imputation.rmdup --exclude merge-merge.missnp --make-bed --out chr22.imputation.rmdup.rm3
plink --bfile chr22.chip.rmdup.flip --exclude merge-merge.missnp --make-bed --out chr22.chip.rmdup.flip.rm3
plink --bfile chr22.imputation.rmdup.rm3 --bmerge chr22.chip.rmdup.flip.rm3 --make-bed --out merge
plink --bfile merge  --genome --out merge.ibd
```
