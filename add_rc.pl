### add readcounts for ref and var to maf file ##

#!/usr/bin/perl
use strict;
use warnings;
(my $usage = <<OUT) =~ s/\t+//g;
This script will add readcounts for ref and var to maf file 
perl run_dir f_maf f_out
OUT

die $usage unless @ARGV == 3;
my ($run_dir,$f_maf,$f_out)=@ARGV;
my $sn_1="";
 
my $sn_2="";
 
my %n_ref=(); 
my %n_var=(); 
my %t_ref=(); 
my %t_var=(); 
my $l; 
open(OUT,">$f_out");

foreach my $l (`cat $f_maf`) 
	{

	my $ltr=$l; 
	chomp($ltr);
	my @t=split("\t",$ltr); 
 
	if($ltr=~/^Hugo/ || $ltr=~/^#version/) { if($ltr=~/^Hugo/) { print OUT $ltr,"\n"; } }

	else {

	 	$sn_1=$t[15]; $sn_1=~s/_T//g; 
	        #print $sn_2,"\t",$sn_1,"\n"; 
		#<STDIN>;		
		if($sn_2 eq "" || $sn_2 ne $sn_1)
		{

		 %n_ref=();
  		 %n_var=();
		 %t_ref=();
 		 %t_var=();
		 $sn_2=$sn_1; 		 

		 my $f_rc_n=$run_dir."/".$sn_1."/".$sn_1.".N.rc.vaf"; 
		 my $f_rc_t=$run_dir."/".$sn_1."/".$sn_1.".T.rc.vaf";
 
		 open(INn,"<$f_rc_n"); 
		 open(INt,"<$f_rc_t");

		 my $id; 

		 my $n_ref;
		 my $n_var; 
		 my $t_ref; 
		 my $t_var;
 
		 while(<INn>)
  		 {
    		   $l=$_;
    		   chomp($l);
		   if($l=~/^#/) { next; } 
		   my @t2=split("\t",$l);
    		   my $id=$t2[1]."_".$t2[2]."_".$t2[4]."_".$t2[5];
		   $n_ref{$id}=$t2[7];
		   $n_var{$id}=$t2[8]; 
		 }	
		
		while(<INt>)
		{
		   $l=$_;
                   chomp($l);
                   if($l=~/^#/) { next; }
                   my @t2=split("\t",$l);
                   my $id=$t2[1]."_".$t2[2]."_".$t2[4]."_".$t2[5];
                   $t_ref{$id}=$t2[7];
                   $t_var{$id}=$t2[8];		  
		}
		}

		my $id2=$t[4]."_".$t[5]."_".$t[10]."_".$t[12]; 	

		#print $id2,"\t",$t_var{$id2},"\t",$n_var{$id2},"\t",$sn_1,"\n"; 
		#<STDIN>;
 	
		if(defined $t_var{$id2} && (length($t[10])<=100) && (length($t[12])<=100))
		{ 

		print OUT $t[0]; 
		for(my $i=1;$i<39;$i++)
		{ 
		print OUT "\t",$t[$i]; 
		}

		if(defined $n_var{$id2})
		{
		print OUT "\t",$t_ref{$id2}+$t_var{$id2},"\t",$t_ref{$id2},"\t",$t_var{$id2},"\t",$n_ref{$id2}+$n_var{$id2},"\t",$n_ref{$id2},"\t",$n_var{$id2}; 						
		}

  		else
		{
		 print OUT "\t",$t_ref{$id2}+$t_var{$id2},"\t",$t_ref{$id2},"\t",$t_var{$id2},"\t","0","\t","0",,"\t","0","\n";
		} 

		for(my $i=45;$i<scalar @t;$i++) 
		{
		print OUT "\t",$t[$i]; 		
		}
		print OUT "\n"; 
		}	
		
		else { 
		  print $id2,"\n";
		  print $sn_1,"\n";	
			}
	}
	}

