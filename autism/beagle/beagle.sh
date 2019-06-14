cd /gpfs/home/guosa/hpc/autism/data/2lof/beagle

mkdir temp
wget https://faculty.washington.edu/browning/beagle/beagle.16May19.351.jar

for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=8 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Djava.io.tmpdir=./temp/ -Xmx32g -jar beagle.16May19.351.jar impute=false gt=All_samples_Exome_QC.temp.vcf.recode.clean.chr$i.vcf.recode.vcf ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz map=~/hpc/db/hg19/beagle/plink.chr$i.GRCh37.map out=All_samples_Exome_QC.chr$i.vcf >>$i.job
qsub $i.job
done
