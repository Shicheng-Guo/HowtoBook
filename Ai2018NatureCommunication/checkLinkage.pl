
use Cwd;
use strict;
my $dir=getcwd;
chdir $dir;
my $threshold=shift @ARGV;
my @snp=glob("rs*.ld");
foreach my $snp(@snp){
open F,$snp;
my ($rs)=split/\./,$snp;
while(<F>){
chomp;
my $line=$_;
my (undef,$chr,$rs1,undef,undef,$rs2,$r,$dp)=split/\s+/,$line;
next if !/$rs/;
if($dp>=$threshold){
        print "$rs\t$line\n";
        last;
}
}
close F;
}
