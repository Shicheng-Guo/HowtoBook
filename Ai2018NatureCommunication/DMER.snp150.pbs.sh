
for i in `ls *bed`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=16 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -q longq  >> $i.job
echo cd /gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication >> ${i}.job
echo bedtools intersect -wao -a $i -b ~/hpc/db/hg19/commonsnp150.hg19.bed  > $i.cSNP150.bed >>  ${i}.job
echo ${i}.job
qsub $i.job
done

