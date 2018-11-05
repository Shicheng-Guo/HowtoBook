use strict;
use Cwd;
use POSIX;
my $dir = getcwd;
chdir "/home/local/MFLDCLIN/guosa/hpc/db/hg19/GTEx_Analysis_v7_eQTL";
my $cpgsnp="/home/local/MFLDCLIN/guosa/hpc/db/hg19/plan2/hg19.CpGSNP.bed4";
my @file =glob("*.v7.egenes.txt");
my %cpgsnp;

open OUT, ">CpG-SNP_eQRL.bed";
open F1,$cpgsnp;
while(<F1>){
        chomp;
        my @line=split/\s+/;
    $cpgsnp{$line[3]}=$_;
}
foreach my $file(@file){
my ($tissue,undef)=split/\.v7/,$file;
open F2,$file;
my $N;
my $J;
while(<F2>){
        chomp;
        next if /gene_id/;
        my @line=split/\s+/;
        $N++;
        next if $line[28]>0.05;
        if($cpgsnp{$line[18]}){
                $J++;
                print OUT "$cpgsnp{$line[18]}\t$line[1]\t$tissue\t$line[15]\t$line[16]\t$line[21]\t$line[28]\n";
                }
        }
print "$J\t$N\t$tissue\n";
}
close OUT;
