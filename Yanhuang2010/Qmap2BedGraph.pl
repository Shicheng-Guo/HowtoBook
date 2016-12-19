#!/usr/bin/perl
use strict;

my $cpgfile=shift @ARGV;
my %cpgpos;
open F,$cpgfile || die "cannot open $cpgfile\n";
while(<F>){
        chomp;
        my($chr,$pos)=split/\s+/;
        $cpgpos{$chr}{$pos}++;
        #print "$chr\t$pos\n";
}
close F;

my @file=glob("*chr*.txt.bedgraph");
my %methdata;
foreach my $file(@file){
        open F,$file || die "cannot open $file\n";
        while(<F>){
        my ($chr,$start,$end,$mf)=split/\s+/;
        $methdata{$chr}{$start}=$mf;
        # print "$chr\t$start\t$end\t$mf\n";
        }
}

foreach my $chr(sort keys %cpgpos){
        foreach my $pos(sort keys %{$cpgpos{$chr}}){
        my $start=$pos-1;
        my $end=$pos;
        next if ! defined $methdata{$chr}{$pos};
        print "$chr\t$start\t$end\t$methdata{$chr}{$pos}\n";
        }
}
