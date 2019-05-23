use strict;
open F,"input.txt";
my %status;
my %input;

while(<F>){
        chomp;
        my ($pid,undef,undef,$status,$num)=split/\"/;
        $pid=substr($pid,0,9);
        $input{$pid}{$status}=$num;
        $status{$status}=$status;
}
my $header=join("\t",sort keys %status);
print "\t$header\n";
foreach my $pid(sort keys %input){
        print "$pid";
        foreach my $status(sort keys %status){
        if(defined $input{$pid}{$status}){
        print "\t$input{$pid}{$status}";
        }else{
        print "\t";
        }
        }
        print "\n";
}
