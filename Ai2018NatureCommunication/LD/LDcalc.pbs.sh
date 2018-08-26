cd /gpfs/home/guosa/nc/plink
awk '{print $1,2*$2-$7,2*$2-$6,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.mirror.hg19.bed 
bedtools intersect -wao -a DMER.mirror.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.mirror.hg19.SNP150.bed
awk '{print $1,$6,$7,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.hg19.bed 
bedtools intersect -wao -a DMER.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.hg19.SNP150.bed
mkdir temp
mkdir LD
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

