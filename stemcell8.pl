use strict;
use LWP::Simple;
use LWP 5.64;
use HTTP::Cookies;
use Encode qw(encode decode);
my $browser=LWP::UserAgent->new();
foreach my $i(1..69){
my $url="http://www.stemcell8.cn/forum-75-$i.html";
my $response=$browser->get($url);
my $content=$response->content;
my @content=split /\n/,$content;
foreach my $tmp(@content){
	if($tmp=~/>(PDF.+)</i){
	print "$1\n";	
	}
}
}
