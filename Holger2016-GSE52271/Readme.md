# Epigenomic analysis detects aberrant super-enhancer DNA methylation in human cancer
## Dataset Analysis: GSE52270,GSE52271 and GSE52272
```bash
cd /home/shg047/oasis/Estellar2016/sortbam
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ../mhb/WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.cor.bed > saminfo.txt
perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt non bisreadMapper /home/shg047/oasis/db/hg19/hg19.chrom.sizes /home/shg047/oasis/db/hg19/HsGenome19.CpG.positions.txt
/home/shg047/oasis/Estellar2016/hapinfo
```
