#!/usr/bin/perl
# USAGE: perl qmap2bedgraph.pl GSE17972_HUMtg5lib.qmap.chr21.txt > GSE17972_HUMtg5lib.qmap.chr21.txt.bedgraph
# Tool: transfer qmap bis-seq format To bedgraph
# Shicheng Guo
# 12/19/2016

use strict;

my $file=shift @ARGV;
my $chr;
if($file=~/(chr\w+)/){
        $chr=$1;
}

open F,$file || die "cannot open $file\n";
my @data=<F>;
foreach my $i(0..$#data){
        my ($read1,$M1,$UM1)=split/\s+/,$data[$i];
        my ($read2,$M2,$UM2)=split/\s+/,$data[$i+1];
        if($read1 eq "C" && $read2 eq "G" ){
                next if (($M1+$M2+$UM1+$UM2)<5);
                my $mf=sprintf("%.3f",($M1+$M2)/($M1+$M2+$UM1+$UM2));
                my $start=$i-1;
                my $end=$start+1;
                print "$chr\t$start\t$end\t$mf\n";
        }
}
