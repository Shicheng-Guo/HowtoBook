open F, "../newphen.txt";
while(<F>){
my($fn1,$fn2)=split/\s+/;
my @file=glob("$fn1*");
foreach my $file(@file){
my $orign=$file;
$file=~ s/$fn1/$fn2/g;
system("cp $orign ./new/$file");
}
}
