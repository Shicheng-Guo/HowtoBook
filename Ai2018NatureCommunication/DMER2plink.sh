rm ./plink/*
for i in chr{1..22} chrX chrY
do
for j in `ls *.pair.bed`
do
grep  -w $i $j >> ./plink/$j.$i.plink.bed
done
done

