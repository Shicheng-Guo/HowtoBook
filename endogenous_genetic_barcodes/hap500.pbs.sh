for i in  ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo sh ./hap500.sh $i 12 135006516 60181 >>$i.job
qsub $i.job
done

# echo sh ./hap500.sh $i 1 249250621 10177 >>$i.job
# echo sh ./hap500.sh $i 2 243199373 10179 >>$i.job
# echo sh ./hap500.sh $i 3 198022430 60069 >>$i.job
# echo sh ./hap500.sh $i 4 191154276 10005 >>$i.job
# echo sh ./hap500.sh $i 5 180915260 10043 >>$i.job
# echo sh ./hap500.sh $i 6 171115067 63854 >>$i.job
# echo sh ./hap500.sh $i 7 159138663 14808 >>$i.job
# echo sh ./hap500.sh $i 8 155270560 11740 >>$i.job
# echo sh ./hap500.sh $i 9 146364022 10163 >>$i.job
# echo sh ./hap500.sh $i 10 141213431 60494 >>$i.job
# echo sh ./hap500.sh $i 11 135534747 61395 >>$i.job
# echo sh ./hap500.sh $i 12 135006516 60181 >>$i.job
# echo sh ./hap500.sh $i 13 133851895 19020047 >>$i.job
# echo sh ./hap500.sh $i 14 115169878 19000017 >>$i.job
# echo sh ./hap500.sh $i 15 107349540 20000041 >>$i.job
# echo sh ./hap500.sh $i 16 102531392 60086 >>$i.job
# echo sh ./hap500.sh $i 17 90354753 52 >>$i.job
# echo sh ./hap500.sh $i 18 81195210 10083 >>$i.job
# echo sh ./hap500.sh $i 19 78077248 60842 >>$i.job
# echo sh ./hap500.sh $i 20 63025520 60343 >>$i.job
# echo sh ./hap500.sh $i 21 59373566 9411239 >>$i.job
# echo sh ./hap500.sh $i 22 59128983 16050075 >>$i.job
# echo sh ./hap500.sh $i X 51304566 60020 >>$i.job
# echo sh ./hap500.sh $i Y 48129895 2655180 >>$i.job
