#  germline_calling_v1.0.pl #
#  the input bam without readgroup information #	
### Song Cao ###
### last updated 7/6/2017 ###
### last updated 10/09/2017 ###
### last updated 10/19/2017 ###
#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;
my $version = 1.0;
#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";
#usage information

(my $usage = <<OUT) =~ s/\t+//g;
This script will process germline callings. 
Pipeline version: $version
$yellow     Usage: perl $0  --srg --step --sre --rdir --ref --log $normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<srg> = bam having read group or not: 1, yes and 0, no (default 1)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<sre> = re-run: 1, yes and 0, no  (default 0)
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
with chr: /gscmnt/gc3027/dinglab/medseq/fasta/GRCh37V1/GRCh37-lite-chr_with_chrM.fa
without chr: /gscmnt/gc3027/dinglab/medseq/fasta/GRCh37/GRCh37-lite.fa
mmy: /gscmnt/gc2737/ding/Reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.20170530.fa 
hg19: /gscmnt/gc2521/dinglab/cptac3/ref/Homo_sapiens_assembly19.fasta

<run_folder> = full path of the folder holding files for this sequence run
<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$red      	 [1]  Run gatk
$red         [2]  Run varscan
$red 		 [3]  Run Pindel
$yellow 	 [4]  Parse pindel
$purple 	 [5]  Merge calls
$green 		 [6]  VCF2MAF
$cyan 		 [7]  Generate final maf
$normal
OUT

#die $usage unless @ARGV == 2;
#my ( $run_dir, $step_number ) = @ARGV;
#if ($run_dir =~/(.+)\/$/) {
#    $run_dir = $1;
#}
#die $usage unless ($step_number >=0)&&(($step_number <= 10));
#GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip
# obtain script path
#my $run_script_path = `dirname $0`;
#__DEFAULT NUMBER OF BINS IE (MUST BE INTEGER)
my $step_number = -1;
my $status_rg = 1;
my $status_rerun=0;
#__HELP (BOOLEAN, DEFAULTS TO NO-HELP)
my $help = 0;

#__FILE NAME (STRING, NO DEFAULT)
my $run_dir="";
my $log_dir="";
my $h37_REF="";

#__PARSE COMMAND LINE
my $status = &GetOptions (
      "step=i" => \$step_number,
      "srg=i" => \$status_rg,
      "sre=i" => \$status_rerun,
      "rdir=s" => \$run_dir,
      "ref=s"  => \$h37_REF,
      "log=s"  => \$log_dir,
      "help" => \$help,
    );

#print $status,"\n";

if ($help || $run_dir eq "" || $log_dir eq ""  || $step_number<0 || $step_number>8) {
      print $usage;
      exit;
   }

print "run dir=",$run_dir,"\n";
print "step num=",$step_number,"\n";
print "status rerun=",$status_rerun,"\n";
print "status readgroup=",$status_rg,"\n";

if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

my $email = "scao\@wustl\.edu";
# everything else below should be automated
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];
my $HOME1=$log_dir;
#store job files here
if (! -d $HOME1."/tmpgermline") {
    `mkdir $HOME1"/tmpgermline"`;
}
my $job_files_dir = $HOME1."/tmpgermline";
#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR_GERMLINE") {
    `mkdir $HOME1"/LSF_DIR_GERMLINE"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR_GERMLINE";
#GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip
# obtain script path
my $script_dir="/gscuser/scao/scripts/git/germlinewrapper";

my $run_script_path=$script_dir;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";
print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
my $h37_REF_bai=$h37_REF.".fai";
my $gatk="/gscuser/scao/tools/GenomeAnalysisTK.jar";
my $STRELKA_DIR="/gscmnt/gc2525/dinglab/rmashl/Software/bin/strelka/1.0.14/bin";
#my $h37_REF="/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37V1/GRCh37-lite-chr_with_chrM.fa";
my $f_exac="/gscmnt/gc2741/ding/qgao/tools/vcf2maf-1.6.11/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
#my $h37_REF_bai="/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37/GRCh37-lite-chr_with_chrM.fa.fai";
my $pindel="/gscuser/qgao/tools/pindel/pindel";
my $PINDEL_DIR="/gscuser/qgao/tools/pindel";
#my $gatk="/gscuser/scao/tools/GenomeAnalysisTK.jar";
my $gatkexe="/gscmnt/gc2525/dinglab/rmashl/Software/bin/gatk/3.7/GenomeAnalysisTK.jar";
my $picardexe="/gscuser/scao/tools/picard.jar";
my $f_centromere="/gscmnt/gc3015/dinglab/medseq/Jiayin_Germline_Project/PCGP/data/pindel-centromere-exclude.bed";
my $java_dir="/gscuser/scao/tools/jre1.8.0_121";

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure
#&check_input_dir($run_dir);
# start data processsing

if ($step_number < 7) {
    #begin to process each sample
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
                $current_job_file="";
                if($step_number==0)
                {
                   &bsub_gatk();
                   &bsub_varscan();
                   &bsub_pindel();
				   &bsub_parse_pindel();
				   &bsub_merge_vcf();
				   &bsub_vcf_2_maf();
                   #&bsub_vep();
                }
                 elsif ($step_number == 1) {
                    &bsub_gatk();
                } elsif ($step_number == 2) {
                    &bsub_varscan(1);
                }elsif ($step_number == 3) {
                    &bsub_pindel(1);
                }elsif ($step_number == 4) {
                    &bsub_parse_pindel(1);
                }elsif ($step_number == 5) {
                    &bsub_merge_vcf(1);
                }elsif ($step_number==6) {
				    &bsub_vcf_2_maf(1);		
				}
				}
				}
		}
	}

if($step_number==7 || $step_number==0) {

# src/generate_report.pm

}


