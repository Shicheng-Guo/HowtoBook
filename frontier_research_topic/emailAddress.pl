#!/usr/bin/perl -w
use strict;
my $SearchResultFile = shift @ARGV;
open FH, $SearchResultFile or die "Cannot open the $SearchResultFile:$!\n";
open OUT, ">$SearchResultFile_Author_Address.xls" or die "Cannot create the output file:$!\n";
print OUT "Email\tTitle\tFirst Name\tMiddle Name\tLast Name\tContribution Title\r\n";
my %AuthorInfo;
my ($tmpLastName, $tmpForeName, $tmpEmail);
while (<FH>) {
  if (m{LastName>(\w.*)</LastName>}){
    $tmpLastName = $1;
  }elsif (m{ForeName>(\w.*)</ForeName}){
    $tmpForeName = $1;
  }elsif (m{<Affiliation>}){
    if (m{\s(\S+?@\w.*)\.</Aff}){
      $tmpEmail = $1;
    }else{
      $tmpEmail = "None";
    }
    if ($tmpEmail =~ m/@/){
      print OUT "$tmpEmail\tDr\t$tmpForeName\t\t$tmpLastName\t\r\n";
      $tmpEmail = "None";
    }
  }
}
