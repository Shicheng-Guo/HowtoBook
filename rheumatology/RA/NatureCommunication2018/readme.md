dss2bedgraph and then upload to ucsc
```
mkdir temp
for i in `ls *dss`
do
echo \#PBS -N $i  > ./temp/$i.job
echo \#PBS -l nodes=1:ppn=1 >> ./temp/$i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> ./temp/$i.job
echo \#PBS -m abe  >> ./temp/$i.job
echo \#PBS -o $(pwd)/temp/ >>./temp/$i.job
echo \#PBS -e $(pwd)/temp/ >> ./temp/$i.job
echo cd $(pwd) >> ./temp/$i.job 
echo awk \'\$3\>0{print \$1,\$2-1,\$2,\$4\/\$3}\' OFS=\"\\t\" $i \> $i.bed >> ./temp/$i.job
qsub ./temp/$i.job
done

bedtools sort -i GSE112658.RA.hg19.bedgraph > GSE112658.RA.sort.hg19.bedgraph &
bedtools sort -i GSE112658.OA.hg19.bedgraph > GSE112658.OA.sort.hg19.bedgraph &

cat GSE112658.RA.sort.hg19.bedgraph >> RA.hg19.bedgraph
cat GSE112658.OA.sort.hg19.bedgraph >> OA.hg19.bedgraph

gzip -c RA-FLS.hg19.bedgraph > RA-FLS.hg19.bedgraph.gz &
gzip -c OA-FLS.hg19.bedgraph > OA-FLS.hg19.bedgraph.gz &

bedtools intersect -wo -a RA-Naturecommnunication.DMR.bed -b ~/hpc/db/hg19/refGeneV2.hg19.bed > RA-Naturecommnunication.DMR.RefGene.hg19.bed
awk '{print $15"."$3,$1,$2,$3,$15"."$3}' OFS="\t" RA-Naturecommnunication.DMR.RefGene.hg19.bed | sort -u > RA-Naturecommnunication.DMR.RefGene.hg19.uni.bed
```
