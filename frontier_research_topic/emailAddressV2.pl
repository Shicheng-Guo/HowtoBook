#!/usr/bin/perl -w
use strict;
my $SearchResultFile = shift @ARGV;
open FH, $SearchResultFile or die "Cannot open the $SearchResultFile:$!\n";
open OUT, ">$SearchResultFile.Author_Address.txt" or die "Cannot create the output file:$!\n";
print OUT "Email\tTitle\tFirst Name\tMiddle Name\tLast Name\tContribution Title\r\n";
my %AuthorInfo;
my %email;
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
	  next if $tmpEmail=~/@\w.*\s*\w.*@/;
	  next if $tmpEmail=~/,/;
	  next if $tmpEmail=~/\s+/;
	  my @email=split/\./,$tmpEmail;
      print OUT "$tmpEmail\tDr\t$tmpForeName\t$tmpLastName\t$tmpLastName\tResearch\t$email[$#email]\n" if ! defined $email{$tmpEmail};
	  $email{$tmpEmail}=$tmpEmail;
      $tmpEmail = "None";
    }
  }
}
