#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;

my $file=shift @ARGV;
open F,$file;
while(<F>){
chomp;
my($chr,$pos,$meth,$unmeth)=split /\t/;
my $start=$pos-2;
my $end=$pos-1;
my $ratio=$meth/($meth+$unmeth+0.001);
print "chr$chr\t$start\t$end\t$ratio\n";
}
