```
cd /gpfs/home/guosa/hpc/project/pmrp/phase2/phase
for i in {1..24}
do
plink --bfile S_Hebbring_Unr.Guo --recode vcf --chr $i --snps-only just-acgt --out ./beagle/S_Hebbring_Unr.Guo.Forward.chr$i
done
```
