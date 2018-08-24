
use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my $genome="/gpfs/home/guosa/hpc/db/hg19/1000Genome";
mkdir "Shuffle" if ! -e "Shuffle";
my $file="/gpfs/home/guosa/nc/GWAS-378.DMER.bed.sort.bed.distance.bed";

open F,$file || die "cannot open $file\n";
my %snp;
while(<F>){
my ($chr,$start1,$end1,$rs,undef,$start2,$end2,undef,undef,$epi)=split/\s+/;
my(undef,$newchr)=split/chr/,$chr;
my($oldfrag,$newfrag,$tmp);
if($start1<=$start2){
$oldfrag="$newchr:$start1-$end2";
$tmp=2*$start1-$start2;
$newfrag="$newchr:$tmp-$start1";
}else{
$oldfrag="$newchr:$start2-$end1";
$tmp=2*$end1-$end2;
$newfrag="$newchr:$start1-$tmp";
}
open OUT,">./Shuffle/$chr.$rs.job";
print OUT "#PBS -N $rs\n";
print OUT "#PBS -o ./Shuffle/temp/\n";
print OUT "#PBS -e ./Shuffle/temp/\n";
print OUT "cd $dir/Shuffle\n";
print OUT "tabix -h $genome/$chr.vcf.gz $oldfrag > $chr.$rs.1.vcf\n";
print OUT "tabix -h $genome/$chr.vcf.gz $newfrag > $chr.$rs.2.vcf\n";
print OUT "plink --vcf $chr.$rs.1.vcf --r inter-chr dprime --out $chr.$rs.1\n";
print OUT "plink --vcf $chr.$rs.2.vcf --r inter-chr dprime --out $chr.$rs.2\n";
close OUT;
system("qsub ./Shuffle/$chr.$rs.job")
}
