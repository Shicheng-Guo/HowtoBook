#!/usr/bin/perl
# USAGE: perl qmap2bedgraph.pl GSE17972_HUMtg5lib.qmap.chrY.txt.gz > GSE17972_HUMtg5lib.qmap.chrY.bedgraph
# Tool: transfer qmap bis-seq format To bigwig/bedgraph/wig
# Shicheng Guo
# Oct/31/2016

use strict;
my $f=shift @ARGV;
my $pos=0;
my $chr;
open F,$f;
if($f=~/(chr\w+)/){
$chr=$1;
}
while(<F>){
chomp;
next if /^\s+$/;
my @line=split/\t/;
my $M=$line[1];
my $UM=$line[2];
$M=0 if !defined $M;
$UM=0 if !defined $UM;
if(($M+$UM)>=0){
my $mf=sprintf("%.3f",$M/($M+$UM+0.1));
my $start=$pos;
my $end=$pos+1;
print "$chr\t$start\t$end\t$M\t$UM\t$mf\n";
}
$pos++;
}

