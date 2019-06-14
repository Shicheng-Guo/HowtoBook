1. Prepare beagle 5.0 and corresponding database 

   https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/APP/beagle5.0/Step1.Prepare.sh

2. Prepare beagle phasing reference for EUR samples (my samples are EUR samples) 

   https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/APP/beagle5.0/Step2.reference.sh

3. Prepare PBS script to run beagle

   https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/APP/beagle5.0/Step3.beaglephasing.sh

4. Post-imputation filiter (DR2>0.6)

   https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/APP/beagle5.0/Step4.postfiltering.sh
   
5. Loss-of-function variants + deleterious SNPs


6. Extract VCF for these LOF Variants 
```
cat All_samples_Exome_QC.chr*.vcf.DR2L0.8.DG.snpEff.High.vcf | grep -v '#' |awk '{print $3}' OFS="\t" > LOF.hg19.txt
cat All_samples_Exome_QC.chr*.vcf.DR2L0.8.DG.vcf.vat.aloft.vcf | grep -v '#' |awk '{print $3}' OFS="\t" | sort -u >>  LOF.hg19.txt
perl -p -i -e 's/chr//' LOF.hg19.txt
sort -u LOF.hg19.txt | grep -v '[.]' > LOF.hg19.sort.txt

for i in {1..22}
do
vcftools --vcf All_samples_Exome_QC.chr$i.vcf.DR2L0.8.vcf --snps LOF.hg19.sort.txt --recode --out samples_Exome_QC.chr$i.vcf.DR2L0.8.LOF.vcf
bcftools annotate -a /gpfs/home/guosa/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') samples_Exome_QC.chr$i.vcf.DR2L0.8.LOF.vcf -o samples_Exome_QC.chr$i.vcf.DR2L0.8.LOF.anno.vcf
done
```

