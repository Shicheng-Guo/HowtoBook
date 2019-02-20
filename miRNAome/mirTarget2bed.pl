use strict;
use Cwd;
use List::Util qw( min max);
chdir getcwd;

open F,shift @ARGV;
while(<F>){
        my ($chr,$start,$end,$gene,$mir)=split/[\s+|:]/;
        $mir=~s/miR-//g;
        if($mir=~/\//){
        my @mir=split/\//g,$mir;
        foreach my $mir(@mir){
        $mir=~s/\.\d//g;
        print "$chr\t$start\t$end\t$gene\thsa-miR-$mir\n";
        }
        }else{
        $mir=~s/\.\d//g;
        print "$chr\t$start\t$end\t$gene\thsa-miR-$mir\n";
        }
}
