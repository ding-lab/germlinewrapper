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

## gatk, vars, pindel ##

my %vaf_rc_gatk=(); 
my %vaf_rc_vs=();
my %vaf_rc_pindel=();


open(OUT,">$f_out");

foreach my $l (`cat $f_maf`) 
	{

	my $ltr=$l; 
	chomp($ltr);
	my @temp=split("\t",$ltr); 
 
	if($ltr=~/^Hugo/) { print OUT $ltr,"\n"; }

	else {

	 	$sn_1=$temp[15]; $sn_1=~s/_T//g; 

		if($sn_2 eq "" || $sn_2 ne $sn_1)
		{

		 %vaf_rc_gatk=();
  		 %vaf_rc_vs=();
		 %vaf_rc_pindel=();
 
		 $sn_2=$sn_1; 		 
		 ##   ## 

		 my $f_gatk_snv=$run_dir."/".$sn_1."/gatk"."/".$sn_1.".snv.gvip.filtered.vcf"; 
		 my $f_gatk_ind=$run_dir."/".$sn_1."/gatk"."/".$sn_1.".indel.gvip.filtered.vcf";
		 my $f_vars_snv=$run_dir."/".$sn_1."/gatk"."/".$sn_1."raw.snp.filtered.vcf";
		 my $f_vars_ind=$run_dir."/".$sn_1."/gatk"."/".$sn_1."raw.indel.filtered.vcf";	
		 my $f_pindel=$run_dir."/".$sn_1."/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf";
 
		 open(INgs,"<$f_gatk_snv"); 
		 open(INgi,"<$f_gatk_ind");
		 open(INvs,"<$f_vars_snv");
		 open(INvi,"<$f_vars_ind");
		 open(INp,"<$f_pindel");

		 my $id; 

		 my $n_ref;
		 my $n_var; 
		 my $t_ref; 
		 my $t_var;
 
		 while(<INgs>)
  		 {
    		my $line=$_;
    		chomp($line);
		    if($line=~/^#/) { next; } 
			my @temp=split("\t",$line);
    		my $id=$temp[0]."_".$temp[1];
    		my $refvar=$temp[3]."_".$temp[4];
    		my $infor=$temp[9];
    		my @temp2=split(":",$infor);
    		my $desc=$temp[8];

    ## remove M allele ##
    		if($refvar=~/M/) { next; }
    		if($refvar=~/R/) { next; }

    		if($desc=~/AD/)
    		{
    		my @temp3=split(",",$temp2[1]);
   	 		my $n_ref=$temp3[0];
    		my $n_var=$temp3[1];
        	$vaf_rc_gatk{$id}=$n_ref."_".$n_var; 
			}
	  	 	
		  }
		  
         close INgs; 
		}	

         while(<INgi>)
         {
            my $line=$_;
            chomp($line);
            if($line=~/^#/) { next; }
            my @temp=split("\t",$line);
            my $id=$temp[0]."_".$temp[1];
            my $refvar=$temp[3]."_".$temp[4];
            my $infor=$temp[9];
            my @temp2=split(":",$infor);
            my $desc=$temp[8];

    ## remove M allele ##
            if($refvar=~/M/) { next; }
            if($refvar=~/R/) { next; }

            if($desc=~/AD/)
            {   
            my @temp3=split(",",$temp2[1]);
            my $n_ref=$temp3[0];
            my $n_var=$temp3[1];
            $vaf_rc_gatk{$id}=$n_ref."_".$n_var;                    
            }

          }
          
         close INgi;
        }

         while(<INvs>)
         {
            my $line=$_;
            chomp($line);

            if($line=~/^#/) { next; }
			
    		my @temp=split("\t",$line);
    		my $id=$temp[0]."_".$temp[1];
    		my $refvar=$temp[3]."_".$temp[4];
    		my $infor=$temp[9];
    		my @temp2=split(":",$infor);
    		my $n_ref_fw=$temp2[-4];
    		my $n_ref_rev=$temp2[-3];
    		my $n_var_fw=$temp2[-2];
    		my $n_var_rev=$temp2[-1];
			my $n_ref_vs=$n_ref_fw+$n_ref_rev; 
			my $n_var_vs=$n_var_fw+$n_var_rev; 
    		if($refvar=~/M/) { next; }
    		if($refvar=~/R/) { next; }
			
			$vaf_rc_vs{$id}=$n_ref_vs."_".$n_var_vs;

		}


         while(<INvi>)
         {
            my $line=$_;
            chomp($line);

            if($line=~/^#/) { next; }

            my @temp=split("\t",$line);
            my $id=$temp[0]."_".$temp[1];
            my $refvar=$temp[3]."_".$temp[4];
            my $infor=$temp[9];
            my @temp2=split(":",$infor);
            my $n_ref_fw=$temp2[-4];
            my $n_ref_rev=$temp2[-3];
            my $n_var_fw=$temp2[-2];
            my $n_var_rev=$temp2[-1];
            my $n_ref_vs=$n_ref_fw+$n_ref_rev;
            my $n_var_vs=$n_var_fw+$n_var_rev;
            if($refvar=~/M/) { next; }
            if($refvar=~/R/) { next; }

            $vaf_rc_vs{$id}=$n_ref_vs."_".$n_var_vs;

        }
### chr1	29902	0/1:20,13 ##
		 while(<INp>)
         {
            my $line=$_;
            chomp($line);

            if($line=~/^#/) { next; }

            my @temp=split("\t",$line);
            my $id=$temp[0]."_".$temp[1];
			
            my $refvar=$temp[3]."_".$temp[4];
            my $infor=$temp[9];
            my @temp2=split(":",$infor);

			my @temp3=split(",",$temp2[1]);
            my $n_ref=$temp3[0];
            my $n_var=$temp3[1];
            $vaf_rc_pindel{$id}=$n_ref."_".$n_var;

            #$vaf_rc_vs{$id}=$n_ref_vs."_".$n_var_vs;

        }		
		
		print OUT $temp[0]; 

		for(my $i=1;$i<39;$i++)
		{ 
		print OUT "\t",$temp[$i]; 
		}

		print OUT "\t",$t_depth,"\t",$t_ref_count,"\t",$t_alt_count,"\t",$n_depth,"\t",$n_ref_count,"\t",$n_alt_count; 						
		for(my $i=45;$i<scalar @temp;$i++) 
		{
		print OUT "\t",$temp[$i]; 		
		}

		print OUT "\n"; 
		
		}	
		
		else { 

		  print $id2,"\n";
		  print $sn_1,"\n";
			
			}
	}
	}
