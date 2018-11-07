
use strict;
use Cwd;
use POSIX;
my $dir = getcwd;
chdir $dir;

open F,"/home/guosa/hpc/db/hg19/allSNP150.hg19.bed";
open OUT,">cpgSNP.hg19.bed.avinput";
while(<F>){
        my ($chr,$start,$end,$rs,$chain,$ref,$alt)=split/\s+/;
        my @allele=split/\//,$alt;
        my (undef,$newchr)=split/chr/,$chr;
        foreach my $allele(@allele){
                next if $allele eq $ref;
                print OUT "$newchr\t$start\t$end\t$ref\t$allele\t$rs\n";
        }
}
close OUT;
