use strict;
use warnings;

my $input=shift @ARGV;
open F,$input;
while(<F>){
	chomp;
	if(/\#\s+ID\s+(.+)/){
	#if(/\#\s+ID\s+(\w+\.\w+\@\w+)/){
	chomp(my $hap1=<F>);
	chomp(my $hap2=<F>);
        my $id=$1;
        $hap1=~s/ //g;
        $hap2=~s/ //g;
	print "$id\t$hap1\t$hap2\n";	
	}
}
