use strict;
use List::Util qw(uniq);
use Cwd;

my @file=glob("out*.uni");
my %data;
my %id;
foreach my $file(@file){
        open F,$file;
        while(<F>){
            chomp;
            my ($gid,$cid,$name,$sex,$age,$dep,$date,$diagnosis,$weight,$height)=split/\t/;
            next if $cid=~/NA/;
            next if length($cid)<18;
			next if $diagnosis=~/\?/;
            $diagnosis=~s/RA|类风湿关节炎/类风湿性关节炎/ig;	
            $diagnosis=~s/OA/骨关节炎/ig;	
			$diagnosis=~s/AS/强直性脊柱炎/ig;	
			$diagnosis=~s/骨质疏松症/骨质疏松/ig;		
            my @disease=split/、|；|;|\s+|,|：|\（|，|\/\\|\(|\//ig,$diagnosis;
            $data{$disease[0]}++;
			$id{$gid}=$gid if ($disease[0] eq "类风湿性关节炎" || $disease[0] eq "骨关节炎" || $disease[0] eq "强直性脊柱炎")
        }
		close F;
}



foreach my $dis(sort keys %data){
next if $data{$dis}<30;
print "$dis\t$data{$dis}\n";
}

foreach my $file(@file){
        open F,$file;
		open OUT,">$file.uni";
        while(<F>){
            my ($gid,undef)=split/\t/;
			print OUT $_ if defined $id{$gid};
		}
		close OUT;
}

a

sub rename{
            my $diagnosis=shift @_;
		    $diagnosis=~s/类风湿关节炎|类风湿性关节炎/类风湿性关节炎/ig;
            $diagnosis=~s/骨关节炎|髌股关节炎|骶髂关节炎/OA/ig;
            $diagnosis=~s/骨质疏松症|骨质疏松/Osteoporosis/ig;
            $diagnosis=~s/系统性红斑狼疮|红斑狼苍/SLE/ig;
            $diagnosis=~s/强直性脊柱炎/AS/ig;	
	        $diagnosis=~s/颈椎骨质增生/CBH/ig;	
            $diagnosis=~s/银屑病关节炎/PSA/ig;	
            $diagnosis=~s/[、?]/;/ig;
			$diagnosis=~s/右膝|左膝|双膝|膝//ig;	
			$diagnosis=~s/早期//ig;	
			$diagnosis=~s/并间质性肺病/-ILD/ig;	
			$diagnosis=~s/甲状腺毒症|甲状腺功能亢进/Thyrotoxicosis/ig;	
			$diagnosis=~s/并间质性肺病/-ILD/ig;	
			$diagnosis=~s/支气管哮喘/BA/ig;	
			$diagnosis=~s/髌|髌股|髌股关节//ig;	
			$diagnosis=~s/髌股关节炎//ig;	
}
