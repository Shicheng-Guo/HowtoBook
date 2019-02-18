##################################################################################
cd /gpfs/home/guosa/hpc/rheumatology/RA/ASA/miRNA-SNP

panel="UTR3miRNAbinding"

awk '{print $3}' UTR.rs.txt | sort -u > 364291.UTRSNP.rsList

for i in {1..22} X Y
do
plink --bfile /gpfs/home/guosa/hpc/db/hg19/1000Genome/chr$i --make-bed --keep ~/hpc/db/hg19/1000Genome/CHB.CHS.txt --allow-no-sex --extract 364291.UTRSNP.rsList --out UTRSNP.CHS.chr$i
plink --bfile UTRSNP.CHS.chr$i --indep-pairphase 500kb 5 0.8 --out UTRSNP.CHS.chr$i
done
cat UTRSNP.CHS.chr*in > UTR3miRNAbinding.ASA.txt
wc -l UTR3miRNAbinding.ASA.txt

for i in {1..22} X Y
do
plink --bfile /gpfs/home/guosa/hpc/db/hg19/1000Genome/chr$i --make-bed --keep ~/hpc/db/hg19/1000Genome/CHB.CHS.txt --allow-no-sex --extract UTR3miRNAbinding.ASA.txt --out UTR3miRNAbinding.ASA.$i
done

awk '{print $1,$4-1,$4,$2}' OFS="\t" UTR3miRNAbinding.ASA.*.bim  >> UTR3miRNAbinding.hg19.bed
wc -l UTR3miRNAbinding.hg19.bed

mkdir temp
## Function Variants in Genome Regions
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo \# bcftools norm -m \+ /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.vcf.bgz -Oz -o gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo \# tabix -p vcf gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -v snps -f PASS -i \'INFO/AF_eas\>0.01 \& INFO/AF_eas\<0.99\' -T $panel.hg19.bed  /gpfs/home/guosa/hpc/db/Gnomad/vcf/gnomad.genomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.genomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.genomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done

ls *genomes*rec.$panel.sort.rmdup.biallelic.vcf.bgz > concat.txt
bcftools concat -f concat.txt -Ov -o gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf
grep -v "#" gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed
awk '{print "chr"$1,$2-1,$2,$3,$4,$5}' OFS="\t" gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.bed > gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.hg19.bed
wc -l gnomad.genomes.r2.1.sites.rec.$panel.merge.vcf.hg19.bed

######################################################################################
