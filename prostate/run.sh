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
