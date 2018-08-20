use Cwd;
use strict;
my $dir=getcwd;
chdir $dir;
my $pop=shift @ARGV;
my $marker=shift @ARGV;
open F,"/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/$marker.bed.sort.bed.SNP150.pair.bed" || die "cannot open $marker\n";
while(<F>){
chomp;
my ($chr,$start,$end,$rs1,$rs2)=split/\s+/;
open OUT,">$marker.$pop.$rs1.$rs2.job";
print OUT "#PBS -N $marker.$pop.$rs1.$rs2\n";
print OUT "#PBS -o ./temp/\n";
print OUT "#PBS -e ./temp/\n";
print OUT "cd $dir\n";
print OUT "plink --vcf /gpfs/home/guosa/nc/$chr.dmer.recode.vcf --keep /gpfs/home/guosa/hpc/db/hg19/1000Genome/$pop.txt --ld $rs1 $rs2 --out ./$pop/$marker.$pop.$rs1.$rs2\n";
print OUT "rm ./$pop/*nosex\n";
close OUT;
system("qsub $marker.$pop.$rs1.$rs2.job");
system("rm $marker.$pop.$rs1.$rs2.job");
}
close F;
