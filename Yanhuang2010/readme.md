#### YanHuang Methylome Project

1. Go to GEO
```
 # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17972
 mkdir /gpfs/home/guosa/hpc/wgbs/GSE17972
 cd  /gpfs/home/guosa/hpc/wgbs/GSE17972
 wget -e robots=off -nH -nd  -r --reject="index.html*" -nd https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17972/suppl/
 rm *index*
 gunzip *.gz
```
2 uncompress and then transfer to bigwig
```
cd /gpfs/home/guosa/hpc/wgbs
for i in {1..22} X Y
do
echo \#PBS -N chr$i  > chr$i.job
echo cd $(pwd) >>chr$i.job
echo gunzip GSE17972_HUMtg5lib.qmap.chr$i.txt.gz >> chr$i.job
qsub chr$i.job
done
```
3. transfer to bigwig
```
for i in `ls *txt`
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo cd /gpfs/home/guosa/hpc/wgbs/GSE17972 >> $i.job
echo ./qmap2bedgraph $i \>$i.bedgraph >> $i.job
echo $i.job
qsub $i.job
done
```
4. txt to bedgraph was suspended since the memory is not big enough, each chrosome requir 30G. I re-do it in CHG1 one by one
```
for i in `ls *txt`
do
perl qmap2bedgraph $i > $i.bedgraph
done

```





##### Li Y, Zhu J, Tian G, Li N et al. The DNA methylome of human peripheral blood mononuclear cells. PLoS Biol 2010 Nov 9;8(11):e1000533. PMID: 21085693
