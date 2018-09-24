use strict;
use Cwd;
chdir getcwd;
my $file="/gpfs/home/guosa/hpc/rheumatology/RA/miRNASNP/hsa.gff3.hg38.txt";
open F,$file || die "cannot open $file";
while(<F>){
        next if /miRNA_primary_transcript/;
        my($chr,undef,undef,$start,$end,undef,$chain,undef,$id)=split/\s+/;
        if($id=~/Name=(.+);/){
        print "$chr\t$start\t$end\t$1\n";
        }
}
close F;
