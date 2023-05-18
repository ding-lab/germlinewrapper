#!/usr/bin/perl

use strict;
use warnings;
die unless @ARGV == 2;

my ($file_in,$dir_out)=@ARGV;

print $file_in,"\n"; 

open(IN,"<$file_in");
#open(OUT,">$file_out"); 
my $f_out; 
my $sn;
my @t; 
my %openfile=();
 
while(<IN>)
{
	my $l=$_; 
	chomp($l);
	if($l=~/^Hugo_Symbol/) { next; } 
	else 
	{ 
		@t=split("\t",$l); 
		$sn=$t[15]; 
		$sn=~s/_T$//g; 
		$f_out=$dir_out."/".$sn."/".$sn.".rc.input.vcf";
		if(!defined $openfile{$sn})
		{
		open(OUT,">$f_out"); 
		$openfile{$sn}=1; 	
		}
		print OUT $t[4],"\t",$t[5],"\t",$t[6],"\t",$t[10],"\t",$t[12],"\n"; 
	}
}

close IN; 
close OUT;
