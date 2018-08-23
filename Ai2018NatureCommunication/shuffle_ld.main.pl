use Cwd;
use strict;
my $dir=getcwd;
chdir $dir;
my $dmer="/gpfs/home/guosa/nc/RA-OA.DMER.GRCH37.bed";
my $gwasnp="/gpfs/home/guosa/nc/GWAS-RA-378.GRCH37.bed";
my $chrsize="~/hpc/db/hg19/hg19.chrom.sizes";
my $wd="/gpfs/home/guosa/nc/ShuffLD";
my $commonSNP150="~/hpc/db/hg19/commonsnp150.hg19.bed";

foreach my $i(1..3){
print "Now wer are on the iteration $i\n";
rmdir "Shuff$i" if ! "Shuff$i";
mkdir "Shuff$i" if ! -e "Shuff$i";
mkdir "./Shuff$i/temp" if ! -e "./Shuff$i/temp";
system("bedtools shuffle -i $dmer -g $chrsize | bedtools sort -i - > ./Shuff$i/dmer.hg19.shuffle.$i.sort.bed");
system("bedtools closest -a $gwasnp -b ./Shuff$i/dmer.hg19.shuffle.$i.sort.bed > ./Shuff$i/dmer.hg19.gwas.shuffle.$i.bed");
&getNearestdmer("./Shuff$i/dmer.hg19.gwas.shuffle.$i.bed","./Shuff$i/dmer.hg19.gwas.shuffle.$i.gwas.bed");
system("bedtools intersect -wao -a ./Shuff$i/dmer.hg19.gwas.shuffle.$i.gwas.bed -b $commonSNP150 > ./Shuff$i/dmer.hg19.gwas.shuffle.$i.gwas.bed.snp150");
open F,"./Shuff$i/dmer.hg19.gwas.shuffle.$i.gwas.bed.snp150" || die "cannot open ";
while(<F>){
        my @line=split/\s+/;
        next if $line[7] eq ".";
        open OUT1,">> ./Shuff$i/$line[0].$line[3].txt";
        print OUT1 "$line[7]\n$line[3]\n";
        close OUT1;
}
close F;

my @snp=glob("./Shuff$i/chr*.txt");
foreach my $snp(@snp){
print "$snp\n";
my @filename=split/\//,$snp;
my ($chr,$rs,undef)=split/\./,$filename[-1];
print "$chr\t$rs\n";
open OUT,">./Shuff$i/$rs.job";
my $pop="~/hpc/db/hg19/1000Genome/ALL2000.txt";
print OUT "#PBS -N $rs\n";
print OUT "#PBS -o ./Shuff$i/temp/\n";
print OUT "#PBS -e ./Shuff$i/temp/\n";
print OUT "cd $dir/Shuff$i/";
print OUT "plink --gzvcf /gpfs/home/guosa/nc/$chr.dmer.recode.vcf.gz --keep $pop --extract $filename[-1] --r inter-chr dprime --out $rs\n";
print OUT "rm *nosex\n";
print OUT "rm *log\n";
close OUT;
system("qsub ./Shuff$i/$rs.job");
}
}

sub getNearestdmer{
my($input,$output) =@_;
open F,$input;
open OUT,">$output";
while(<F>){
        my @line=split/\s+/;
        print OUT "$line[4]\t$line[5]\t$line[6]\t$line[3]\n";
}
close OUT;
}
