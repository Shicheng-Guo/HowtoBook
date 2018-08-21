use Cwd;
use strict;
my $dir=getcwd;
chdir $dir;
my @snp=glob("*.snp.list.uni");
foreach my $snp(@snp){
my ($chr,$rs,undef)=split/\./,$snp;
open OUT,">$rs.job";
print OUT "#PBS -N $rs\n";
print OUT "#PBS -o ./temp/\n";
print OUT "#PBS -e ./temp/\n";
print OUT "cd $dir\n";
print OUT "plink --vcf /gpfs/home/guosa/nc/$chr.dmer.recode.vcf --extract $snp --r inter-chr dprime --out $rs\n";
print OUT "rm *nosex\n";
print OUT "rm *log\n";
close OUT;
system("qsub $rs.job");
}
close F;
