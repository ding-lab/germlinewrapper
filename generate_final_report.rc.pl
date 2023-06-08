#!/usr/bin/perl
use strict;
use warnings;
die unless @ARGV == 1;
my ($run_dir)=@ARGV;

my $working_name= (split(/\//,$run_dir))[-1];

my $f_sum=$run_dir."/".$working_name.".rc.maf\n";
my $f_status=$run_dir."/".$working_name.".rc.status\n";

open(OUT1,">$f_sum");
open(OUT2,">$f_status"); 

my $head_w=0;

foreach my $d (`ls $run_dir`)
{
  	my $dtr=$d; 
	chomp($dtr); 
	my $f_maf=$run_dir."/".$dtr."/".$dtr.".rc.maf"; 
	if(-e $f_maf) 
	{
		my $count=0;
		foreach my $l (`cat $f_maf`) 
		{ 
			my $ltr=$l; 
			chomp($ltr); 
			if(($ltr=~/version/ || $ltr=~/^Hugo/) && $head_w==1)  { next; } 
			else { 
			if($ltr=~/^Hugo/) { $head_w=1; print OUT1 $ltr,"\n"; } 
			else { 
				my @temp=split("\t",$ltr); 
				my $annot=$temp[8];
				if(1==1) 
				       {
					print OUT1 $ltr,"\n"; 
					$count++;
					} 
				} 
        		}
      	}
		print OUT2 $count,"\t",$f_maf,"\n";
	} 
  }


close OUT1;
close OUT2; 
