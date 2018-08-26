use strict;
use Cwd;
use List::Util;

my $dir=getcwd;

chdir $dir;
my (%rsq,%dp);
my @file=glob("*.log");
foreach my $file(@file){
        my($rs1,$rs2,$id,undef)=split/\./,$file;
        open F,$file;
        while(<F>){
                if(/R-sq\s+=\s+(\d+\.\d+)\s+D'\s+=\s+(\d+\.\d+)/){
                push @{$rsq{$rs1}{$id}}=$1;
                push @{$dp{$rs1}{$id}}=$2;
                }
        }
}

open OUT1,">DMER.twoside.R2.txt";
open OUT2,">DMER.twoside.dp.txt";
foreach my $rs(sort keys %rsq){
        print OUT1 "$rs\tmax(@{$rsq{$rs}{"r1"}})\tmax(@{$rsq{$rs}{"r2"}}}\n";
        print OUT2 "$rs\tmax(@{dp{$rs}{"r1"}})\tmax(@{dp{$rs}{"r2"}}}\n";
}
close OUT1;
close OUT2;




