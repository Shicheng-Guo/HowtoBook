#### 2019-06-04

001A_Normal.Mut1        
001A_Normal.Mut2        
001A_Normal.VarD

bgzip 001A_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
bcftools view -s 001A_Normal.Mut1 -o 001A_Normal.Mut1.vcf
bcftools annotate -x INFO,^FORMAT/GT -s 001A_Normal.Mut1 -o 001A_Normal.Mut1.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' 001A_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf.gz

cd /home/local/MFLDCLIN/guosa/hpc/project/LungBrainMetastasis/bed
for i in A B C E F G H I J K L M O
do
wc -l 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed > $i.txt
wc -l 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed >> $i.txt
bedtools intersect -wa -a 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed -b 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed | wc -l >> $i.txt
done

use strict;
my @file=glob("*.txt");
foreach my $file(@file){
open F,$file;
print "$file";
while(<F>){
my ($num)=split/\s+/;
print "\t$num";
}
print "\n"
}

use strict;
my @file=glob("*.bed");
my %cat;
my %type;
foreach my $file(@file){
open F,$file;
my($sam)=split/\_/,$file;
while(<F>){
my @line=split/\s+/;
$cat{$sam}{$line[6]}++;
$type{$line[6]}=$line[6];
}
}

foreach my $sam(sort keys %cat){
print "$sam";
 foreach my $type(sort keys %type){
 print "\t$cat{$sam}{$type}";
 }
 print "\n";
}


cd /home/local/MFLDCLIN/guosa/hpc/project/LungBrainMetastasis/bed
for i in A B C E F G H I J K L M O
do
bedtools intersect -wa -a 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed -b 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.hg19_multianno.txt.bed >$i.overlap.bed
done


#### 2019-07-04



mv S001B_Tumor2-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf   001B_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
mv S001D_Tumor2_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf  001D_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
mv S001L_Tumor2-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001L_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf
mv S001N_Tumor2_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001N_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf

mv S001B_Tumor1-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf    001B_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf
mv S001D_Tumor1_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf  001D_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf
mv S001L_Tumor1-Normal_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001L_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf
mv S001N_Tumor1_mem.merged.HighConf.snpEff_ann.hg19_multianno.vcf 001N_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf

ls -l *Tumor2*.vcf | wc -l
ls -l *Tumor1*.vcf | wc -l

cd /home/guosa/hpc/project/LungBrainMetastasis
awk '{print $1,$2,$4,$5}' OFS="\t" *Normal-Tumor1_mem*.bed  | sort -u > Normal-Tumor1.uni.hg19.bed
awk '{print $1,$2,$4,$5}' OFS="\t" *Normal-Tumor2_mem*.bed  | sort -u > Normal-Tumor2.uni.hg19.bed


## step 1.1 check sample names and remove or keeps samples of interest. 
rm Tumor1_VarD.txt
cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
bcftools query -l 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf | grep Tumor1 | grep VarD >> Tumor1_VarD.txt
done

rm Tumor2_VarD.txt
cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
bcftools query -l 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf | grep Tumor2 | grep VarD >> Tumor2_VarD.txt
done

## Step 2. copy raw vcf to new folder to keep raw data safe
cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
cp 001$i\_Normal-Tumor1_mem.merged.HighConf.snpEff.vcf $i.T1.vcf
done

cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
cp 001$i\_Normal-Tumor2_mem.merged.HighConf.snpEff.vcf $i.T2.vcf
done

## step 3. 
for i in A B C E F G H I J K L M O
do
echo $i
bcftools view --force-samples -S Tumor1_VarD.txt $i.T1.vcf | bcftools annotate -x ID,^INFO/AN,INFO/DP,FORMAT -I +'%CHROM:%POS' | bcftools sort -Ov -o $i.T1.sort.vcf
perl -p -i -e 's/\.\/\./0\/1/g' $i.T1.sort.vcf
bgzip -f $i.T1.sort.vcf
tabix -f -p vcf $i.T1.sort.vcf.gz
done
ls *.T1.sort.vcf.gz > T1.txt
bcftools merge -l T1.txt -Ov | bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') > T1.vcf


cd /home/guosa/hpc/project/LungBrainMetastasis/vcf
for i in A B C E F G H I J K L M O
do
echo $i
bcftools view --force-samples -S Tumor2_VarD.txt $i.T2.vcf | bcftools annotate -x ID,^INFO/AN,INFO/DP,FORMAT -I +'%CHROM:%POS' | bcftools sort -Ov -o $i.T2.sort.vcf
perl -p -i -e 's/\.\/\./0\/1/g' $i.T2.sort.vcf
bgzip -f $i.T2.sort.vcf
tabix -f -p vcf $i.T2.sort.vcf.gz
done
ls *.T2.sort.vcf.gz > T2.txt
bcftools merge -l T2.txt -Ov | bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') > T2.vcf

touch -m T*.vcf

## step 4. Rename chrosome and re-annoate VCF with ANNOVAR

rm rename-chrs.txt
for i in {1..24} X Y
do
echo -e chr$i'\t'$i >> rename-chrs.txt 
done 

bcftools annotate --rename-chrs rename-chrs.txt T1.vcf -Ov -o T1.chr.vcf
bcftools annotate --rename-chrs rename-chrs.txt T2.vcf -Ov -o T2.chr.vcf
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') T1.chr.vcf 
bcftools annotate -a ~/hpc/db/hg19/refGene.hg19.VCF.sort.bed.gz -c CHROM,FROM,TO,GENE -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') T2.chr.vcf 

table_annovar.pl -vcfinput T1.chr.vcf ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 -out T1 -remove -protocol refGene,dbnsfp33a -operation gx,f -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt
table_annovar.pl -vcfinput T2.chr.vcf ~/hpc/tools/annovar/humandb/ --thread 4 -buildver hg19 -out T2 -remove -protocol refGene,dbnsfp33a -operation gx,f -nastring . -otherinfo -polish -xref ~/hpc/tools/annovar/humandb/gene_fullxref.txt
 




