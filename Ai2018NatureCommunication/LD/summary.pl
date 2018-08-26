
use strict;
use Cwd;
my $dir=getcwd;
chdir $dir;
my (%rsq,%dp);
my @file=glob("*.log");
foreach my $file(@file){
        my($rs1,$rs2,$id,undef)=split/\./,$file;
        open F,$file;
        while(<F>){
                if(/R-sq\s+=\s+(\d+\.\d+)\s+D'\s+=\s+(\d+\.\d+)/){
                $rsq{$rs1}{$id}=$1;
                $dp{$rs1}{$id}=$2;
#                print "$rs1\t$id\t$1\t$2\n";
                }
        }
}

open OUT1,">DMER.twoside.R2.txt";
open OUT2,">DMER.twoside.dp.txt";
foreach my $rs(sort keys %rsq){
print OUT1 "$rs\t$rsq{$rs}{'r1'}\t$rsq{$rs}{'r2'}\n";
print OUT2 "$rs\t$dp{$rs}{'r1'}\t$dp{$rs}{'r2'}\n";
}
close OUT1;
close OUT2;
