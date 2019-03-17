sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

open F,"BLM-VEP.hg19.txt";
my %data;
while(<F>){
next if $_ =~ /Location/i;
my(undef,$var,undef,$anno)=split/\s+/;
my @anno=split/\,/,$anno;
push(@{$data{$var}},@anno);
}
foreach my $var(sort keys %data){
        my @anno = uniq @{$data{$var}};
    my $anno=join(";",@anno);
    my ($chr,$pos,undef)=split/[:|-]/,$var;
    my $start=$pos-1;
    my $end=$pos;
    print "chr$chr\t$start\t$end\t$anno\n";
}
