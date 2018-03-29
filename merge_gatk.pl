#!/usr/bin/perl

### snv: union call from varscan and gatk
### indel: callings from pindel or both gatak and varscan

use strict;
use warnings;
die unless @ARGV == 2;

my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
my ($sample_full_path,$sample_name)=@ARGV;
my $f_snv;
my $f_indel; 

my $f_snv_idx; 
my $f_indel_idx; 

my $f_gvip; 
my $f_gvip_idx; 

my $f_snv_out=$sample_full_path."/gatk/".$sample_name.".snv.gvip.vcf";
my $f_ind_out=$sample_full_path."/gatk/".$sample_name.".indel.gvip.vcf";

open(OUT1,">$f_snv_out"); 
open(OUT2,">$f_ind_out"); 

foreach my $chr (@chrlist)
    {

	$f_snv=$sample_full_path."/gatk/".$sample_name.".snv.gvip.$chr.vcf";
	$f_indel=$sample_full_path."/gatk/".$sample_name.".indel.gvip.$chr.vcf";
 
	$f_snv_idx=$sample_full_path."/gatk/".$sample_name.".snv.gvip.$chr.vcf.idx";
	$f_indel_idx=$sample_full_path."/gatk/".$sample_name.".indel.gvip.$chr.vcf.idx";

	$f_gvip=$sample_full_path."/gatk/".$sample_name.".gvip.$chr.vcf"; 
	$f_gvip_idx=$sample_full_path."/gatk/".$sample_name.".gvip.$chr.vcf.idx";
	#$f_raw_idx=$sample_full_path."/gatk/".$sample_name.".gvip.$chr.vcf.idx";	
	#$f_raw=$sample_full_path."/gatk/".$sample_name.".snv.gvip.$chr.vcf";

	if(-e $f_snv) 
	{
		foreach my $l (`cat $f_snv`) 	
		{
		my $ltr=$l; 
		chomp($ltr); 
		if($ltr=~/^#/ &&  !($chr eq "1"))  { next; }
		else { print OUT1 $ltr,"\n"; }
		}
 		`rm $f_snv`;	
		`rm $f_snv_idx`; 
	}

	if(-e $f_indel) 
	{
		foreach my $l (`cat $f_indel`)
    	{
        my $ltr=$l; 
        chomp($ltr);
        if($ltr=~/^#/ &&  !($chr eq "1"))  { next; }
        else { print OUT2 $ltr,"\n"; }
    	}   	
		`rm $f_indel`; 
		`rm $f_indel_idx`; 
 	}

	`rm $f_gvip`; 
	#`rm $f_raw`; 
	#`rm $f_raw_idx`;
	`rm $f_gvip_idx`; 

    #print GATK "rawvcf=".$sample_full_path."/gatk/".$sample_name.".raw.$chr.vcf\n";
    #print GATK "gvipvcf=".$sample_full_path."/gatk/".$sample_name.".gvip.$chr.vcf\n";
    #print GATK "snvvcf=".$sample_full_path."/gatk/".$sample_name.".snv.gvip.$chr.vcf\n";
    #print GATK "indelvcf=".$sample_full_path."/gatk/".$sample_name.".indel.gvip.$chr.vcf\n";
    #print GATK "     ".$run_script_path."genomevip_label.pl GATK \${rawvcf} \${gvipvcf}"."\n";
	}
