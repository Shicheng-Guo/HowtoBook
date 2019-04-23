for i in `ls *.strelka.somatic.vcf.gz`; do bcftools view -f PASS $i | wc -l ; done

for i in `ls *.strelka.somatic.vcf.gz`
do
bcftools view -H -f PASS $i | awk '{print $1,$2,$3,$4,$5}' OFS="\t" >> prostate.vep.input
echo $i
done

for i in `ls *.strelka.somatic.vcf.gz`
do
bcftools view -H -f PASS $i | awk '{print $1,$2-1,$2,$1":"$2-1"-"$2}' OFS="\t" > $i.anno.bed
echo $i
done

perl vep_merge2.pl > prostate.somaticMutation.gene.txt
perl -p -i -e 's/;-//g' prostate.somaticMutation.gene.txt

mkdir temp
for i in `ls *strelka.somatic.vcf.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx4g -jar ~/hpc/tools/snpEff/snpEff.jar "GRCh37.75" $i \> $i.anno.vcf >> $i.job 
qsub $i.job
done






python configureStrelkaSomaticWorkflow.py --normalBam baihaixing_control.recal.bam --tumorBam baihaixing_T1.recal.bam --referenceFasta /thinker/storage/udata/bing/GATK/hg19/ucsc.hg19.fasta --exome --callRegions /thinker/storage/udata/bing/bed/xGenExomeResearchPaneltargets.sorted.ip100.merged.bed.gz --runDir baihaixing_T1_strelka_somatic


mkdir temp
perl -p -i -e "{s/chr//}" $panel.hg19.bed
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
echo \#bcftools norm -m \+ ~/hpc/project/pmrp/phase1/imputation/chr$i.dose.filter.vcf.gz -Oz -o chr$i.rec.vcf.bgz >> $i.job
echo \#tabix -p vcf chr$i.rec.vcf.bgz >> $i.job
echo bcftools view -f PASS -i \'INFO\/STATUS \' ~/hpc/db/Gnomad/vcf/gnomad.exomes.r2.1.sites.chr$i.rec.vcf.bgz -Ou -o  gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz >>$i.job
echo bcftools sort gnomad.exomes.r2.1.sites.chr$i.rec.$panel.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz >> $i.job
echo bcftools norm -d all gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.vcf.bgz -Ou -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz >> $i.job
echo bcftools view -m2 -M2 -v snps gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.vcf.bgz -Ov -o gnomad.exomes.r2.1.sites.chr$i.rec.$panel.sort.rmdup.biallelic.vcf.bgz >>$i.job
qsub $i.job
done
bcftools view -i 'INFO/STATUS ~ "Somatic" | INFO/STATUS ~ "LOH"'  
bcftools view -f ., PASS P2_T3.strelka.somatic.vcf.gz
for i in `ls *.strelka.somatic.vcf.gz`; do bcftools view -f PASS $i | wc -l ; done
for i in `ls *.strelka.somatic.vcf.gz`
do
bcftools view -H -f PASS $i | awk '{print $1,$2,$3,$4,$5}' OFS="\t" >> prostate.vep.input
echo $i
done

### pick up all the somatic mutations from strelka
for i in `ls *.strelka.somatic.vcf.gz`
do
bcftools view -f PASS $i -o $i.pass.vcf 
bgzip $i.pass.vcf
tabix -p vcf $i.pass.vcf.gz
bcftools view -H -f PASS $i | awk '{print $1,$2-1,$2,$1":"$2-1"-",$3,$4,$5}' OFS="\t" > $i.bed
echo $i
done 


### pick up all the somatic mutations from vardict
for i in `ls *.vardict.vcf.gz`
do
bcftools view -i 'INFO/STATUS !="Germline"' $i -o $i.pass.vcf 
bgzip $i.pass.vcf
tabix -p vcf $i.pass.vcf.gz
bcftools view -H -f PASS $i.pass.vcf.gz |awk '{print $1,$2-1,$2,$1":"$2-1"-"$2,$3,$4,$5}' OFS="\t" > $i.bed
echo $i
done 

###### merge vardict and strelka
ls *.vcf.gz.bed | awk -F"." '{print $1}' | sort -u > FileUnit.txt

cat *vardict.strelka.overlap.bed | awk '{print $1,$3}' OFS="\t" | sort -u > valid.pos.txt

for i in `ls *.strelka.somatic.vcf.gz.pass.vcf.gz`
do
bcftools view -T valid.pos.txt $i -Oz -o $i.vardit.overlap.vcf.gz
tabix -p $i.vardit.overlap.vcf.gz
echo $i
done

for i in `ls *.strelka.somatic.vcf.gz.pass.vcf.gz`
do
echo $i
bcftools view -T valid.pos.txt $i -Oz -o $i.vardit.overlap.vcf.gz
tabix -p vcf $i.vardit.overlap.vcf.gz
done

mkdir temp
for i in `ls *strelka.somatic.vcf.gz.pass.vcf.gz.vardit.overlap.vcf.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx16g -jar ~/hpc/tools/snpEff/snpEff.jar "GRCh37.75" $i \> $i.anno.vcf >> $i.job 
qsub $i.job
done

for i in `ls *.strelka.somatic.vcf.gz.pass.vcf.gz.vardit.overlap.vcf.gz.anno.vcf`
do
java -jar ~/hpc/tools/snpEff/SnpSift.jar extractFields $i CHROM POS REF ALT "ANN[*].GENE:" | awk -F'[\t|]' '{print $1,$2,$3,$4,$8}' OFS="\t" > $i.hg19.anno.bed
java -jar ~/hpc/tools/snpEff/SnpSift.jar tstv $i > $i.hg19.tstv.bed
done

for i in `ls *.strelka.somatic.vcf.gz.pass.vcf.gz.vardit.overlap.vcf.gz.anno.vcf`
do
grep -v '#' $i|  perl -lane '{print "@F[0]\t@F[1]\t@F[3]\t@F[4]\t$1" if $_=~/SGT=(\w+->\w+);/}' > $i.hg19.SGT.bed
done

 
for i in P10_T1 P11_T1 P12_T1 P12_T2 P12_T3 P13_T1 P13_T2 P13_T4 P14ming_T1 P14_T1 P14_T3 P14_T4 P15_T1 P16_T1 P16_T2 P1_T1 P1_T2 P2_T2 P2_T3 P2_T4 P3_T2 P3_T3 P3_T4 P4_T1 P4_T3 P5_T1 P5_T2 P5_T3 P6_T1 P6_T2 P6_T4 P7_T1 P8_T1 P8_T2 P9_T2 P9_T3 P9_T4 
do
bcftools view -T $i.vardict.strelka.overlap.bed $i.vardict.vcf.gz.pass.vcf.gz -Oz -o $i.vardict.strelka.overlap.vcf.gz
tabix -p vcf $i.vardict.strelka.overlap.vcf.gz
zcat $i.vardict.strelka.overlap.vcf.gz | grep -v '#' | perl -lane '{print "@F[0]\t@F[1]\t@F[3]\t@F[4]\t$1" if $_=~/TYPE=(\w+);/}' > $i.vardit.strelka.overlap.TYPE.hg19.bed
echo $i
done

mkdir temp
for i in `ls *vardict.strelka.overlap.vcf.gz`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx16g -jar ~/hpc/tools/snpEff/snpEff.jar "GRCh37.75" $i \> $i.anno.vcf >> $i.job 
qsub $i.job
done

for i in `ls *vardict.strelka.overlap.vcf.gz.anno.vcf` 
do
grep -v '#' $i | perl -lane '{print "@F[0]\t@F[1]\t@F[3]\t@F[4]\t$1\t$2\t$3\t$4\t$5" if $_=~/TYPE=(\w+);.*\|(\w+)\|(\w+)\|(\w+)\|(ENSG\w+)/}' > $i.hg19.bed
echo $i
done



