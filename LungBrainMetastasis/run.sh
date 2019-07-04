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




