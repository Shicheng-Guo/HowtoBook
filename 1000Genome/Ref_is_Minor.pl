# https://www.biostars.org/p/119420/

use strict;
use Cwd;
chdir getcwd;
my $file=shift @ARGV;
open F,$file || die "cannot open $file\n";
while(<F>){
next if /^#/;
if(/AF=(\d+\.\d+)/){
my $aaf=$1;
if($aaf>0.5){
my ($chr,$pos,$rs,$ref,$alt)=split/\s+/;
my $start=$pos-1;
print "$chr\t$start\t$pos\t$rs\t$ref\t$alt\t$aaf\n";
}
}
}
