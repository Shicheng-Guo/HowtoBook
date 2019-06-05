### gnomad.genomes.r2.1.sites

cd /home/guosa/hpc/db/Gnomad/vcf
panel="LOF"
mkdir temp
## Function Variants in Exom Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \~ \"stop\" \| INFO\/vep \~ \"lost\" \| INFO\/vep \~ \"splice\" \| INFO\/vep \~ \"gain\" \| INFO\/vep \~ \"frame\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Ov -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf >>$i.job
echo bcftools view -v snps -f PASS -i \'\(INFO\/vep \! ~ \"splice_region\"\)\' ~/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf -Ov -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.final.vcf >>$i.job
qsub $i.job
done
ls *rrec.$panel.final.vcf > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf |grep rs | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.genomes.r2.1.sites.rec.$panel.hg19.vcf.bed
wc -l gnomad.genomes.r2.1.sites.rec.$panel.hg19.vcf.bed

### gnomad.exomes.r2.1.sites

