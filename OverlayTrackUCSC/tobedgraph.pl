#!/usr/bin/perl -w

# Tobedgraph.pl: transfer to bedgraph.(GSE63183)
# Position file is always 1-based.
# Contact: Shicheng Guo
# Version 1.3
# Data: 12/15/2016

use strict;
use warnings;
use Cwd;

use strict;
my %tissue;
my $file=shift @ARGV;

open F2,$file;
my ($sample,undef)=split /.txt|-corrected/,$file;
open OUT,">$file.bedgraph";
print OUT "track type=bedGraph name=\"$sample\" visibility=full color=20,150,20 altColor=150,20,20 windowingFunction=mean\n";
while(<F2>){
chomp;
my ($chr,$pos,$vg,$mf)=split /\t/;
my ($mc,$nmc)=split/\//,$vg;
next if ($mc+$nmc)<5;
my $start=$pos-1;
my $end=$pos;
$mf=sprintf("%.3f",$mf);
print OUT "$chr\t$start\t$end\t$mf\n";
}
