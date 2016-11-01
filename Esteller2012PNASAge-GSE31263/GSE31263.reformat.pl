#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;

my $file=shift @ARGV;
open F,$file;
while(<F>){
chomp;
my($chr,$pos,$pos2,$mf)=split /\t/;
my $start=$pos-1;
my $end=$pos;
my $ratio=$mf;
print "$chr\t$start\t$end\t$ratio\n";
}
