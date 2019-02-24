
for i in `ls *dss`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job
echo awk \'\$3\>4{print \$1,\$2-1,\$2,\$4\/\$3}\' OFS\=\"\\t\" $i \> $i.bedgraph >> ./temp/$i.job
qsub ./temp/$i.job
done

for i in `ls *.bedgraph`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job
echo bedtools sort -i $i \> $i.sort >> ./temp/$i.job
echo bedGraphToBigWig $i.sort /gpfs/home/guosa/hpc/db/hg19/hg19.chrom.sizes $i.bw >> ./temp/$i.job
qsub ./temp/$i.job
done

