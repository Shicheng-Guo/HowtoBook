
#!/usr/bin/perl
use strict;
use Cwd;
chdir getcwd;
open F,shift @ARGV;
while(<F>){
next if /C\/G/;
next if /\/.\//;
next if /-\//;
next if /\/-/;
next if /lengthTooLong/i;
next if /\/\w\w/;
next if /\w\w\//;
print "$_";
}
