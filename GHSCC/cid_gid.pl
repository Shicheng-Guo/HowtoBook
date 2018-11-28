
use strict;
use Cwd;

my @file=glob("out*.txt");
my %data;
foreach my $file(@file){
        open F,$file;
        while(<F>){
                chomp;
            my ($gid,$cid,$name,$sex,$age,$dep,$date,$diagnosis,$weight,$height)=split/\t/;
                push @{$data{$cid}},$gid;
        }
}
foreach my $cid (sort keys %data){
        my $length=length($cid);
        next if ($length < 15);
        next if $cid=~/\./;
        print "$cid";
        foreach my $gid(@{$data{$cid}}){
        print "\t$gid";
        }
        print "\n";
}
