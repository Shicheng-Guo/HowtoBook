#frameshift_variant (not include  inframe_insertion, inframe_deletion )
cd ~/hpc/db/Gnomad/vcf
panel="ExomeFrame"
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"frameshift\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -T ./temp/ >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Oz -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
ls *rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.exomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.exomes.r2.1.sites.rec.$panel.hg19.vcf.bed
