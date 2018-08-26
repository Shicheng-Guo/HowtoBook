cd /gpfs/home/guosa/nc/plink
awk '{print $1,2*$2-$7,2*$2-$6,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.mirror.hg19.bed 
bedtools intersect -wao -a DMER.mirror.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.mirror.hg19.SNP150.bed
awk '{print $1,$6,$7,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.hg19.bed 
bedtools intersect -wao -a DMER.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.hg19.SNP150.bed
mkdir temp
mkdir LD
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed


# Whole Sample 
cd /gpfs/home/guosa/nc/plink
awk '{print $1,2*$2-$7,2*$2-$6,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.mirror.hg19.bed 
bedtools intersect -wao -a DMER.mirror.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.mirror.hg19.SNP150.bed
awk '{print $1,$6,$7,$4}' OFS="\t" ../GWAS-378.DMER.bed.sort.bed.distance.bed >  DMER.hg19.bed 
bedtools intersect -wao -a DMER.hg19.bed -b ~/hpc/db/hg19/commonsnp150.hg19.bed > DMER.hg19.SNP150.bed
mkdir temp
mkdir LD
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --ld "$4" "$8 " --out './ld/'"$4"."$8".r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

 
# CHB Sample
cd /gpfs/home/guosa/nc/plink
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CHB.r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CHB.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CHB.r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

# CEU Sample
cd /gpfs/home/guosa/nc/plink
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CEU.r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/CEU.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".CEU.r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

# YRI Sample
cd /gpfs/home/guosa/nc/plink
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/YRI.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".YRI.r1 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.mirror.hg19.SNP150.bed 
awk '$8!="." {cmd="plink --bfile ~/hpc/db/hg19/1000Genome/"$1 " --keep '/gpfs/home/guosa/hpc/db/hg19/1000Genome/YRI.txt' --ld "$4" "$8 " --out './ld/'"$4"."$8".YRI.r2 | qsub -N "$4"."$8 " -e ./temp/ -o ./temp/";system(cmd)}' DMER.hg19.SNP150.bed

