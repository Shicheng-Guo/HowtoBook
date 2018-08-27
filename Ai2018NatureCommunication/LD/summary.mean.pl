
use strict;
use Cwd;
use List::Util qw/max sum/;
use Statistics::Basic qw(:all);
my $dir=getcwd;
chdir $dir;
my (%rsq,%dp);
my @file=glob("*.log");
foreach my $file(@file){
        my($rs1,$rs2,$id,undef)=split/\./,$file;
        open F,$file;
        while(<F>){
                if(/R-sq\s+=\s+(\d+\.\d+)\s+D'\s+=\s+(\d+\.\d+)/){
                push @{$rsq{$rs1}{$id}},$1 if $1;
                push @{$dp{$rs1}{$id}},$2 if $2;
                print "$rs1\t$rs2\t$id\t$1\t$2\n";
                }
        }
}
open OUT1,">DMER.twoside.R2.mean.txt";
open OUT2,">DMER.twoside.dp.mean.txt";
foreach my $rs(sort keys %rsq){
        my $r1max=max @{$rsq{$rs}{"r1"}};
        my $r2max=max @{$rsq{$rs}{"r2"}};
        my $dp1max=mean @{$dp{$rs}{"r1"}};
        my $dp2max=mean @{$dp{$rs}{"r2"}};
        print OUT1 "$rs\t$r1max\t$r2max\n" if defined $r1max && defined $r2max;
        print OUT2 "$rs\t$dp1max\t$dp2max\n" if defined $dp1max && defined $dp2max;
}
close OUT1;
close OUT2;
