#!/usr/bin/perl
use strict;
my $input=shift @ARGV;
open F,$input;
print "\{\n";
print "\t\"ids\":[\n";
while(<F>){
next if /filename/;
my($uuid,undef)=split/\s+/;
print "\t\"$uuid\",\n";
}
print "\t\]\n";
print "\}"
