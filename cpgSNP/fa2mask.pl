use strict;
use Cwd;
use POSIX;
my $dir = getcwd;
chdir $dir;

my %iupac=(
'-/C' => 'Q',
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
my $input="~/hpc/db/hg38/commonSNP150.hg38";
my $chr=shift @ARGV;
sub match_all_positions {
	my @regex1= qw /YG MG SG HG BG VG NG QG/;
	my @regex2= qw /CR CK CS CB CV CD CN/;
    my @regex3;
	for my $i(qw /Y M S H B V N Q/){
		for my $j(qw /R K S B V D N/){
			push @regex3,"$i$j";
		}
	}
	my $regex1=join "|",@regex1;
	my $regex2=join "|",@regex2;
	my $regex3=join "|",@regex3;
    my ($string) = shift @_;
    my @ret;
    while ($string =~ /$regex1/g) {
        push @ret, [$-[0], $+[0],"C"];
    }
    while ($string =~ /$regex2/g) {
        push @ret, [$-[0], $+[0],"G"];
    }
    while ($string =~ /$regex3/g) {
        push @ret, [$-[0], $+[0],"C"];
    }
    return @ret
}

sub match_CpG_positions {
	my @regex= qw /CG/;
	my $regex=join "|",@regex;
    my ($string) = shift @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, [$-[0], $+[0],"CG"];
    }
    return @ret
}

open F1,"$chr.fa";
my $genome;
while(<F1>){
    chomp;
    next if />/;
    $_=~s/N/X/g;
    $genome .=$_;
}
my @seq=split //,$genome;

open F2,"$input" || die "cannot open $input! (commonSNP or allSNP)\n";
open OUT,">$chr.mask.fa";
my %database;
while(<F2>){
    chomp;
    my $line=$_;
    my (undef,$chr,$start,$end,$rs,undef,$strand,undef,$ref,$obs)=split/\s+/,$line;
    my $position=$end-1;
	next if !$iupac{$obs};
    $seq[$position]=$iupac{$obs};
    $database{"$end"}=$line;
}
close F2;

my $cline=0;
foreach my $seq(@seq){
	$cline++;
	print OUT "$seq";
	print OUT "\n" if $cline % 100 ==0;
}
close OUT;
