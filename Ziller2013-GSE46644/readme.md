## Ziller MJ, Gu H, Müller F, Donaghey J et al. Charting a dynamic DNA methylation landscape of the human genome. Nature 2013 Aug 22;500(7463):477-81. PMID: 23925113


```bash
cd /oasis/tscc/scratch/shg047/Ziller2013/fastq
perl ~/bin/smartbismark.pl --input saminfo.txt --genome hg19 --server TSCC --submit no --queue hotel
for i in SRR949210 SRR949211 SRR949212 SRR949213 SRR949214 SRR949215
do
qsub $i.pbs
done
```
