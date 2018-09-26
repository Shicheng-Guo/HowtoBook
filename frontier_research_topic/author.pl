open F,"xml.xml";
while(<F>){

        if(/<AuthorList/){
        chomp(my $f1=<F>);
                chomp(my $f2=<F>);
                chomp(my $f3=<F>);
                chomp(my $f4=<F>);
                chomp(my $f5=<F>);
                chomp(my $f6=<F>);
        if($f3=~/\>(.+)\</){
                print "$1\t";
        }
        if($f2=~/>(.+)</){
                print "$1\t";
        }
        if($f6=~/(\W+\@\w+.\w+)/){
                print "$1";
        }
        print "\n";
        }
}
