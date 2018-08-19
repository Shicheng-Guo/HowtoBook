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

```






##### Li Y, Zhu J, Tian G, Li N et al. The DNA methylome of human peripheral blood mononuclear cells. PLoS Biol 2010 Nov 9;8(11):e1000533. PMID: 21085693
