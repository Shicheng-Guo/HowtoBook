#!/usr/bin/perl
# USAGE: perl qmap2bedgraph.pl GSE17972_HUMtg5lib.qmap.chrY.txt.gz > GSE17972_HUMtg5lib.qmap.chrY.bedgraph
# Tool: transfer qmap bis-seq format To bigwig/bedgraph/wig
# Shicheng Guo
# Oct/31/2016

use strict;
my $f=shift @ARGV;
my $pos;
my $chr;
open F,$f;
if($f=~/(chr\w+)/){
$chr=$1;
}
while(<F>){
chomp;
my @line=split/\t/;
my $M=$line[1];
my $UM=$line[2];
if(($M+$UM)>=5){
my $mf=sprintf("%.3f",$M/($M+$UM));
my $start=$pos;
my $end=$pos+1;
print "$chr\t$start\t$end\t$mf\n";
}
$pos++;
}
