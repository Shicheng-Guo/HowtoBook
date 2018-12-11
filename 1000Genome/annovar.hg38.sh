wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi
gunzip ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

# split to each chrosome to speed up annotation
cd /gpfs/home/guosa/hpc/db/hg38/1000genome
for i in {1..22} X Y
do
echo \#PBS -N $i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo vcftools --gzvcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --chr $i --recode --out ./annovar/chr22 >>chr$i.job
qsub chr$i.job
done

# download reference database for annnovar hg38
cd ~/hpc/tools/annovar/
annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/ &
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/ &
annotate_variation.pl -downdb -webfrom annovar exac03 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2_all humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar gnomad_exome humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar 1000g2015aug humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar kaviar_20150923 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar hrcr1 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cg69 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar gnomad_genome humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar dbnsfp30a humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar esp6500siv2 humandb -buildver hg38 &
annotate_variation.pl -downdb esp6500siv2 humandb -buildver hg38 &
annotate_variation.pl -downdb -webfrom annovar gerp++ humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cadd humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cadd13 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar fathmm humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar eigen humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar gwava humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar dbscsnv11 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar spidex humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar clinvar_20160302 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar cosmic70 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar icgc21 humandb -buildver hg38  &
annotate_variation.pl -downdb -webfrom annovar nci60 humandb -buildver hg38  &


# annotation with ANNOVAR from vcf to excel for each chrosome.
cd /gpfs/home/guosa/hpc/db/hg38/1000genome
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo \#PBS -l nodes=1:ppn=1 >> chr$i.job
echo cd $(pwd) >> chr$i.job
echo vcftools --gzvcf ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --chr $i --recode --out ./annovar/chr$i >>chr$i.job
echo convert2annovar.pl -format vcf4 -allsample -withfreq -includeinfo chr$i.recode.vcf  \> ../annovar/chr$i.vcf.avinput >> chr$i.job
echo table_annovar.pl ../annovar/chr$i.vcf.avinput ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg38 --csvout -out ../annovar/chr$i.hg38 -remove -protocol refGene,dbnsfp35c,gwasCatalog,wgRna -operation gx,f,r,r -nastring . -otherinfo  -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt >> chr$i.job
qsub chr$i.job
done
rm *.e*
rm *.o*
