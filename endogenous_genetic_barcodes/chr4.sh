cd /gpfs/home/guosa/hpc/db/hg19/1000Genome

for i in `seq 10300 250 191154276`
do
j=$((i+500))
tabix -p vcf ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 4:$i-$j > ./hap500/temp.$i.vcf
cut -f 10- ./hap500/temp.$i.vcf > ./hap500/temp.$i.txt
perl -p -i -e "s/\|/\t/g" ./hap500/temp.$i.txt
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' ./hap500/temp.$i.txt > ./hap500/temp.$i.trans.txt
sort -u ./hap500/temp.$i.trans.txt | wc -l >> chr4.hap.txt
done
