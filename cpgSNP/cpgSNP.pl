use strict;
use Cwd;
use POSIX;
my $dir = getcwd;
chdir $dir;

chomp(my $chr=shift @ARGV);

my %iupac1=(
'A/G' => 'R',
'C/T' => 'Y',
'A/C' => 'M',
'G/T' => 'K',
'C/G' => 'S',
'A/T' => 'W',
'A/C/T' => 'H',
'C/G/T' => 'B',
'A/C/G' => 'V',
'A/G/T' => 'D',
'A/C/G/T' => 'N',
);

my %iupac2=(
'A/G' => 'Y',
'C/T' => 'R',
'A/C' => 'K',
'G/T' => 'M',
'C/G' => 'S',
'A/T' => 'W',
'A/C/T' => 'D', 
'C/G/T' => 'V',
'A/C/G' => 'B', 
'A/G/T' => 'H',
'A/C/G/T' => 'N',
);

sub match_all_positions {
	my @regex1= qw /YG MG SG HG BG VG NG/;
	my @regex2= qw /CR CK CS CB CV CD CN/;
    my @regex3;
	for my $i(qw /Y M S H B V N/){
		for my $j(qw /R K S B V D N/){
			push @regex3,"$i$j";
		}
	}
	my $regex1=join "|",@regex1;
	my $regex2=join "|",@regex2;
	my $regex3=join "|",@regex3;
    my ($string) = shift @_;
    my @ret;
    while ($string =~ /$regex1/ig) {
        push @ret, [$-[0], $+[0],"C"];
    }
    while ($string =~ /$regex2/ig) {
        push @ret, [$-[0], $+[0],"G"];
    }
    while ($string =~ /$regex3/ig) {
        push @ret, [$-[0], $+[0],"C"];
    }
    return @ret
}

open F1,"$chr.fa"|| die "cannot open $chr.fa!\n";
my $genome;
while(<F1>){
    chomp;
    next if />/;
    $_=~ s/N/X/ig;
    $genome .=$_;
}
my @seq=split //,$genome;
close F1;

open F2,"/gpfs/home/guosa/hpc/db/hg38/commonSNP150.hg38" || die "cannot open ~/db/hg38/commonSNP150.hg38\n";
my %database;
while(<F2>){
	next if /bin/;
    chomp;
    my $line=$_;
    my (undef,$mychr,$start,$end,$rs,undef,$strand,undef,$ref,$obs,undef)=split/\s+/,$line;
	next if $mychr ne $chr;
        my $position=$end-1;
	if($strand eq "+"){
	next if !$iupac1{$obs};
	print "$mychr\t$chr\t$position\t$rs\t$ref\t$obs\t$iupac1{$obs}\n";
        $seq[$position]=$iupac1{$obs};
        $database{$end}=$line;
	}else{
	next if !$iupac2{$obs};
	print "$mychr\t$chr\t$position\t$rs\t$ref\t$obs\t$iupac2{$obs}\n";
        $seq[$position]=$iupac2{$obs};
        $database{$end}=$line;
	}
}
close F2;

open OUT,">$chr.mask.fa";
my $cline=0;
foreach my $seq(@seq){
	$cline++;
	print OUT "$seq";
	print OUT "\n" if $cline % 200 ==0;
}
close OUT;

my $maskgenome=join "",@seq;
open OUT2,">$chr.CpGSnp.bed" || die "cannot open $chr.CpGSnp.bed!\n";
my @pos=&match_all_positions($maskgenome);
foreach my $pos(@pos){
	my $start= $$pos[2] eq "C"? $$pos[0]:($$pos[0]+1);
	my $end =$$pos[2] eq "C"? $$pos[0]+1:($$pos[0]+2);
	my $mode=join "", @seq[($$pos[0]-3)..(($$pos[1]+1))];
	my @output=split /\s+/,$database{$end};
	my $LossGain= $output[8] eq "C" ? "Loss" : "Gain";
	print OUT2 "$chr\t$start\t$end\t$output[4]\t$output[6]\t$LossGain\t$output[8]\t$output[9]\t$mode\n";
}
close OUT2;
