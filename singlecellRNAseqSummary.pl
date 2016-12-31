# This is for Song to summary his single-cell RNA-seq dropseq data
# 12/16/2015

#!/usr/bin/perl
my $file=shift @ARGV;
open F,$file;
my %counts;
while(<F>){
chomp;
my (@A)=split/\t/;
foreach my $a(@A[3..9]){
        next if ! defined $a;
        $counts{$A[0]}{$a}++;
}
}

print "\t1\t2\t3\t4\t5\t6\t7\t8\n";
foreach my $key(sort keys %counts){
open OUT,">summary.txt";
print "$key";
foreach my $i(qw/1 2 3 4 5 6 7 8/){
my $value=$counts{$key}{$i};
$value="NA" if ! defined($value);
print "\t$value"
}
print "\n";
}
