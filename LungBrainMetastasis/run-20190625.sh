/home/guosa/hpc/project/LungBrainMetastasis/bcftools

001O_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
bcftools annotate -x ID, INFO, FORMAT LungBrain.vcf.gz -Oz -o LungBrain.trim.vcf.gz
bcftools annotate -x ID,FORMAT,INFO,FILTER 001O_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf | grep 29486505


 
bcftools annotate -Ov -I +'%ID' #leaves it as the existing ID
bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' #sets it to chr:pos:ref:alt



for i in `ls *snpEff.vcf`
do
bcftools norm -m-any $i | bcftools norm -Ov --check-ref w -f ~/hpc/db/hg19/hg19.fa | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS:%REF:%ALT' -Oz -o $i.trim.vcf.gz
done

120714401
242383170
105181012
7977928
77216330
70317174
45664573
8093805
45664573
29486505
26534837

chr17:80396771

for i in `ls *snpEff.vcf`
do
bcftools view -T ^remove.txt $i -Ov | bcftools norm -m-any -Ov --check-ref w -f ~/hpc/db/hg19/hg19.fa | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS' -Oz -o $i.trim.vcf.gz
done

for i in `ls *.gz`
do
bcftools query -l $i | grep VarD | grep Tumor
done

for i in `ls *snpEff.vcf`
do
bcftools view --force-samples -T ^remove.txt -S sample.txt $i -Ov | bcftools norm -m-any -Ov --check-ref w -f ~/hpc/db/hg19/hg19.fa | bcftools annotate -x ID,INFO,FORMAT -I +'%CHROM:%POS' -Oz -o $i.trim.vcf.gz
done

for i in `ls *.gz`
do
done

for i in `ls *.merged.HighConf.snpEff.vcf.trim.vcf.gz`
do
bcftools view --force-samples -S sample.txt $i -Ov -o $i.sin.vcf
done

for i in `ls *.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf`
do
perl updata.pl $i > $i.chan.vcf
done

for i in `ls *.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf.chan.vcf`
do
bcftools sort $i -Oz -o $i.sort.vcf.gz
done

for i in `ls *mem.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf.chan.vcf.sort.vcf.gz`
do
tabix -f -p vcf $i
done


bcftools query -l LungBrainAnnovar.hg19_multianno.trim.vcf | grep Tumor1 > Lung.txt
bcftools query -l LungBrainAnnovar.hg19_multianno.trim.vcf | grep Tumor2 > Brain.txt

bcftools view -S Lung.txt LungBrainAnnovar.hg19_multianno.vcf -Ov -o LungBrainAnnovar.Lung.hg19.vcf
bcftools view -S Brain.txt LungBrainAnnovar.hg19_multianno.vcf -Ov -o LungBrainAnnovar.Brain.hg19.vcf

bcftools annotate -x "^INFO/AN,INFO/AC,INFO/GENE,INFO/ExonicFunc.refGene" LungBrainAnnovar.Lung.hg19.vcf | less -S 

for i in `ls *.vcf`
do 
grep -v '#' $i | wc -l 
done

for i in `ls *mem.merged.HighConf.snpEff.vcf.trim.vcf.gz.sin.vcf.chan.vcf`
do 
echo $i
# grep -v '#' $i | wc -l 
done
