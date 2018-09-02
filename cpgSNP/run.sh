cd /home/guosa/hpc/db/hg19/plan2/commonCpGSNP
rm cpgSNP.hg19.uni.bed
for i in `ls *cpgsnp.bed`
do
perl unitrim.pl $i >> cpgSNP.hg19.uni.bed
done
awk '{print $4}' cpgSNP.hg19.uni.bed > cpgSNP.hg19.list.txt 

cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
plink --bfile /home/guosa/hpc/db/hg19/1000Genome/chr$i --extract /home/guosa/hpc/db/hg19/plan2/commonCpGSNP/cpgSNP.hg19.list.txt --make-bed --out chr$i.cpgSNP
done

cd /home/guosa/hpc/cpgSNP
for i in {1..22} X Y
do
perl assignRefence.pl chr$i.cpgSNP.bim > chr$i.cpgSNP.bim.tr
mv chr$i.cpgSNP.bim.tr chr$i.cpgSNP.bim
awk '{print $2"\tC"}' chr$i.cpgSNP.bim > chr$i.C.ref
plink --bfile chr$i.cpgSNP --reference-allele cpgSNP.ref.txt --make-bed --out chr$i.cpgSNP.C
done
plink --bfile /home/guosa/hpc/cpgSNP/chr5.cpgSNP.C --chr 5 --from-bp 317989 --to-bp 325025 --recode fastphase --reference-allele cpgSNP.ref.txt  --freq --out /home/guosa/hpc/cpgSNP/AHRR

