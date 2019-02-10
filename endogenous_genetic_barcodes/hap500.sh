#!/bin/bash 
for i in `seq $4 250 $3`
do
j=$((i+500))
tabix -h -p vcf ALL.chr$2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $2:$i-$j > ./hap500/temp.chr$2.$1.$i.vcf
bcftools view -H ./hap500/temp.chr$2.$1.$i.vcf -S $1.1000G.S1.txt > ./hap500/temp.chr$2.$1.$i.vcf.txt
rm ./hap500/temp.chr$2.$1.$i.vcf
cut -f 10- ./hap500/temp.chr$2.$1.$i.vcf.txt > ./hap500/temp.chr$2.$1.$i.txt
rm ./hap500/temp.chr$2.$1.$i.vcf.txt
perl -p -i -e "s/\|/\t/g" ./hap500/temp.chr$2.$1.$i.txt
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
}' ./hap500/temp.chr$2.$1.$i.txt > ./hap500/temp.chr$2.$1.$i.trans.txt
rm ./hap500/temp.chr$2.$1.$i.txt
sort -u ./hap500/temp.chr$2.$1.$i.trans.txt | wc -l >> chr$2.$1.$1.hap.txt
rm ./hap500/temp.chr$2.$1.$i.trans.txt
done

