rm rename-chrs.txt
for i in {1..24} X Y
do
echo -e chr$i'\t'$i >> rename-chrs.txt 
done 

bcftools annotate --rename-chrs rename-chrs.txt T1.vcf
