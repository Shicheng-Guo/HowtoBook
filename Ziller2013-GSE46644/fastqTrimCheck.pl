#!/usr/bin/perl -w

# Script to check which samples are failed in Fastq Trim with trim_galore
# Contact: Shihcheng.Guo@Gmail.com
# Version 1.3
# input: redundant bedgraph
# output: non-redundant bedgraph
# Run the script in the fold of coverage files created by bismark alignmetor.

use strict;
use warnings;
use Cwd;
use File::Basename;
use Term::ANSIColor;
 
my $dir=getcwd;
chdir "../fastq";
my @fastq=glob("*.fastq.gz");
chdir "../fastq_trim";
my @fastq_trim=glob("*.fq.gz");

my @fastq2=map{SampleName($_)} map{basename($_,".gz")}@fastq;
my @fastq_trim2=map{SampleName($_)} map{basename($_,".fq.gz")}@fastq_trim;

foreach my $i (0..$#fastq2){
$fastq_trim2[$i]="NA" if ! defined $fastq_trim2[$i];
print "$fastq2[$i]\t$fastq_trim2[$i]\n";
}

compare(\@fastq2,\@fastq_trim2);

sub SampleName{
my($SampleName,undef)=split/.fastq|_val_/,$_;	
return($SampleName)
}

sub compare{
my($Fastq,$FastqTrim)=@_;
foreach my $fastq (@$Fastq){
unless( $fastq ~~ @$FastqTrim){
print color("red"), "$fastq missed in fastq_trim folder\n",color("reset");
}
}
}

