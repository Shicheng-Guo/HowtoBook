use strict;
use Cwd;
chdir "/gpfs/home/guosa/hpc/db/Gnomad";
open F,shift @ARGV;
while(<F>){
next if /^#/;
if(/vep=T\|(\w+)/){
my @line=split/\s+/;
print "$line[0]\t$line[1]\t$line[1]\t$line[3]\t$line[4]\t$1\t$line[2]\n";
}
}
