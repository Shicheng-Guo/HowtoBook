#!/usr/bin/perl -w

# A perl script to merge bismark cov file by SRX list.
# Contact: Shihcheng.Guo@Gmail.com
# Version 1.3
# Go to http://sra.dnanexus.com/studies/SRP028600/samples
# Select SRS and Click Related RUNS then get the Table as the input
# Run the script in the fold of coverage files created by bismark alignmetor.

use strict;
use warnings;
use Cwd;

my $file=shift @ARGV;
my %SRA;
open F,$file;
while(<F>){
chomp;
if(/(SRR\d+)/){
	my $SRR=$1;
	if(/(SRX\d+)/){
		my $SRX=$1;
		print "$SRR\t$SRX\n";
		push @{$SRA{$SRX}},$SRR;
		}
	}
}

my @cov=glob("*.cov");

my %mf;
my %pos;

foreach my $SRX(sort keys %SRA){
        foreach my $SRR (@{$SRA{$SRX}}){
        	foreach my $cov(@cov){
        		if($cov=~/$SRR/){
        			open F,$cov;
        			while(<F>){
        				chomp;
        				my($chr,$start,$end,undef,$nme,$nume)=split/\s+/;
        				my $pos="$chr:$start-$end";
 						#print "$pos\t$nme\t$nume\n";
                        $pos{$SRX}{$pos}=$pos;
                        $mf{$SRX}{$pos}{'me'} += $nme;
                        $mf{$SRX}{$pos}{'ume'}+= $nume;
        			}
        		}
        	}
        }
}

foreach my $SRX(sort keys %SRA){
	open OUT,">$SRX.bedgraph";
	foreach my $pos(sort keys %{$pos{$SRX}}){
		my $nme=$mf{$SRX}{$pos}{'me'};
        my $nume=$mf{$SRX}{$pos}{'ume'};
        if(($nme+$nume)>=5){
        my ($chr,$start,$end)=split/[:-]/,$pos;
        my $mf=$nme/($nme+$nume+0.01);
        print OUT "$chr\t$start\t$end\t$mf"; 	
        }
	}
}







