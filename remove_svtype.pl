#!/usr/bin/perl

use strict;
use warnings;
die unless @ARGV == 2;
my ($file_in,$file_out)=@ARGV;

open(IN,"<$file_in");
open(OUT,">$file_out"); 

while(<IN>)
{
	my $l=$_; 
	chomp($l);
	if($l=~/^#/) { print $l, "\n"; next; }

	$l=~s/SVTYPE=//g; 
	print OUT $l,"\n"; 
}

close IN; 
close OUT;
