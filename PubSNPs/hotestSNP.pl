#/home/guosa/hpc/Genotype/SNP.txt
use strict;
use Cwd;
chdir getcwd;
my $file=shift @ARGV;
my %snp;
open F,$file || die "cannot open $file\n";
open OUT,">SNP.count.txt";
my $j;
while(<F>){
if(/(rs\d+)/){
$snp{$1}++;
}
}
foreach my $snp(sort keys %snp){
print OUT "$snp\t$snp{$snp}\n";
}
close F;
close OUT;
