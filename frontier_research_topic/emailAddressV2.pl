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
	  next if $tmpEmail=~/\>/;
	  next if $tmpEmail=~/Affiliatio/;
	  next if $tmpEmail=~/\s+/;
	  next if $tmpEmail=~/mail\:/;
	  next if $tmpEmail=~/\:/;
	  next if $tmpEmail=~/\(/;
	  next if $tmpEmail=~/\)/;
	  next if $tmpEmail=~/\//;
	  next if $tmpEmail=~/\\/;
	  next if $tmpEmail=~/gov1/;
	  next if $tmpEmail=~/stanfordedu/;
      $tmpEmail=lc $tmpEmail;
	  my @email=split/\./,$tmpEmail;
	  next if defined $email{$tmpEmail};
      print OUT "$tmpEmail\tDr\t$tmpForeName\t$tmpLastName\t$tmpLastName\tResearch\n" if ! defined $email{$tmpEmail};
	  $email{$tmpEmail}++;
      $tmpEmail = "None";
	  my $country=$email[$#email];
    }
  }
}

foreach my $emailtemp(sort keys %email){
print "$emailtemp\t$email{$emailtemp}\n";
}
