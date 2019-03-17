sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

open F,"BLM-VEP.hg19.txt";
my %data;
while(<F>){
my(undef,$var,undef,$anno)=split/\s+/;
my @anno=split/\,/,$anno;
push(@{$data{$var}},@anno);
}
foreach my $var(sort keys %data){
        my @anno = uniq @{$data{$var}};
    my $anno=join(";",@anno);
    my ($chr,$start,$end)=split/[:|-]/,$var;
    print "chr$chr\t$start\t$end\t$anno\n";
}
