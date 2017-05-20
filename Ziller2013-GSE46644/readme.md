## Charting a dynamic DNA methylation landscape of the human genome

```bash
cd /oasis/tscc/scratch/shg047/Ziller2013/fastq
perl ~/bin/smartbismark.pl --input saminfo.txt --genome hg19 --server TSCC --submit no --queue hotel
# Primary Colon Cancer: SRR949210 SRR949211 SRR949212 
# Adjacent normal tissue: SRR949213 SRR949214 SRR949215
for i in SRR949210 SRR949211 SRR949212 SRR949213 SRR949214 SRR949215
do
qsub $i.pbs
done
```

### Description

Paired tumor and normal colon whole-genome bisulfite sequencing (WGBS) data. More information about this data is available on GEO ([GSE46644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46644)). 

### Citation

Ziller MJ, Gu H, Müller F, Donaghey J et al. Charting a dynamic DNA methylation landscape of the human genome. Nature 2013 Aug 22;500(7463):477-81. PMID: 23925113

