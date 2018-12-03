
PTPN11: rs3750050, rs9640663
SNP1="rs3750050"
SNP2="rs9640663"
Genesymbol="PTPN11"
grep ~/hpc/db/hg19/PTPN11
tabix -fh ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 7:67256713-87256713 > $Genesymbol.vcf
plink --vcf $Genesymbol.vcf --maf 0.01 --make-bed --keep ~/hpc/db/hg19/1000Genome/CHS.txt --freq --allow-no-sex --out $Genesymbol
grep $SNP1 $Genesymbol.frq
grep $SNP2 $Genesymbol.frq
plink --bfile $Genesymbol --ld $SNP1 $SNP2
plink --bfile $Genesymbol --ld $SNP1 $SNP2 
plink --vcf $Genesymbol.vcf --maf 0.01 --make-bed --keep ~/hpc/db/hg19/1000Genome/CHB.txt --freq --out $Genesymbol
grep $SNP1 $Genesymbol.frq
grep $SNP2 $Genesymbol.frq
