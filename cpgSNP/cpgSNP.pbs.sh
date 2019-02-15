#####################################################################
## Function Variants in Genome Regions
# /gpfs/home/guosa/hpc/db/hg19/
wget https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/cpgSNP/cpgSNP.pl
mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=16 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >>./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job
echo perl cpgSNP.pl chr$i >> ./temp/$i.job
qsub ./temp/$i.job
done
#########################################################################
