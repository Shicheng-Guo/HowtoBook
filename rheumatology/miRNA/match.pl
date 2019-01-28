use strict;
use Cwd;
use POSIX;
my $dir = getcwd;
chdir $dir;

open F1,"miRNA.UTR.rs.sort.txt";
my %data;
while(<F1>){
	chomp;
	my ($gene,$mirna,$rs,$type)=split/\s+/;	
    $data{"$gene#$mirna"}{$type}=$rs;
}

foreach my $key(sort keys %data){
	my ($gene,$mirna)=split/\#/,$key;
	print "$gene\t$mirna\t$data{$key}{'miRNA'}\t$data{$key}{'UTR'}\n" if defined $data{$key}{'miRNA'} && $data{$key}{'UTR'};
}
