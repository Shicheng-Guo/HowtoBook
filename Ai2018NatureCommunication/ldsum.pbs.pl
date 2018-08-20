use Cwd;
use strict;
my %db;
my $dir=getcwd;
chdir $dir;
my $marker=shift @ARGV; # H3K4me1.bed.sort.bed.SNP150.pair.bed
my $pop=shift @ARGV;
chdir "./$pop";
my @file=glob("$marker.*.*.*.log");
foreach my $file(@file){
open F,$file;
my($marker,$pop,$rs1,$rs2)=split/\./,$file;
while(<F>){
chomp;
if(/D\'\s+=\s+(\d\.\d+)/){
$db{$pop}{$marker}{$rs1}{$rs2}=$1;
$db{$pop}{$marker}{$rs2}{$rs1}=$1;
}
}
close F;
}

open F,"/gpfs/home/guosa/hpc/rheumatology/RA/NatureCommunication/$marker.bed.sort.bed.SNP150.pair.bed" || die "cannot open $marker\n";
open OUT,">../$marker.$pop.dmer.ld.ld.txt";
while(<F>){
my $r2;
chomp;
my($chr,$start,$end,$rs1,$rs2)=split/\s+/;
if($db{$pop}{$marker}{$rs1}{$rs2}){
        $r2= $db{$pop}{$marker}{$rs1}{$rs2};
}else{
        $r2="NA";
}
print OUT "$_\t$r2\n";
}
close F;
close OUT;
