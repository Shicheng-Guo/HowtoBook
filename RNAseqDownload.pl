#!/usr/bin/perl

# Download RNA-seq data from SRA and then rename the name for STAR analysis

use strict;
use Cwd;
chdir getcwd;
open F,"/oasis/tscc/scratch/zhl002/RNAseq/Kolodziejczyk_data/E-MTAB-2600.sdrf.txt";
while(<F>){
next if /Characteristics/;
my @line=split/\t/;
my $ftp=$line[26];
my @ftp=split /[\/|\_|\.]/,$ftp;
print "$ftp\t$ftp[11]\t$ftp[12]\t$ftp[13]\n";
my $sampleid=$ftp[11];
my $reads=$ftp[13];
if(! -d $ftp){
mkdir $sampleid;
}
if($reads eq "1"){
system("wget $ftp -O /home/shg047/oasis/Alice/RNAseq/$sampleid/$sampleid\_L001_R1_001.fastq.gz");
}elsif($reads eq "2"){
system("wget $ftp -O /home/shg047/oasis/Alice/RNAseq/$sampleid/$sampleid\_L001_R2_001.fastq.gz");
}
}
