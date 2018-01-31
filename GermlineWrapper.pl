#  GermlineWrapper.pl
#  the input bam without readgroup information #	
### Song Cao and Matthew Wyczalkowski ###

# To do: implement command line arguments:
# -d: dry run.  print commands without executing them

#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;
my $version = 2.0;
#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";
#usage information

require('src/run_GATK.pl');
require('src/run_varscan.pl');
require('src/run_pindel.pl');
require('src/parse_pindel.pl');
require('src/merge_vcf.pl');

# TODO: parameters should be passed in configuration file
# TODO: test for existence of BAM, reference, etc.

(my $usage = <<OUT) =~ s/\t+//g;
This script will perform germline calls for WXS and WGS data
$yellow     Usage: perl $0 step_number config_file [config_file_2] $normal
step_number run this pipeline step by step. (running the whole pipeline if step number is 0)
config_file Input configuration file.  See below for format
config_file_2 Optional secondary configuration file, any parameters here override configuration previous configuration

$red      	 [1]  Run gatk
$red         [2]  Run varscan
$red 		 [3]  Run Pindel
$yellow 	 [4]  Parse pindel
$purple 	 [5]  Merge calls
$normal

Config File details

Format:
    key = value

** TODO **
* Incorporate original germline wrapper arguments
* confirm all optional config file parameters are used
* confirm that status_rg is necessary and implemented correctly. 
    This is a pre-processing step and may be better to remove it
* do we need path to sw_dir?  Or to germline wrapper installation?  May need both
* find/replace all "somatic" instances

Required configuration file keys
    normal_bam
    reference_fasta
    sample_name

Optional configuration file parameters  - NOTE, these are taken from SomaticWrapper, must be cleaned up
    assembly - GRCh37 or GRCh38.  Used only for VEP, not currently implemented
    sw_dir - Somatic Wrapper installation directory.  Default is /usr/local/somaticwrapper  *** TODO change this ***
    sw_data - Somatic Wrapper analysis results directory.  Default is /data/data
            Per-sample analysis directory is sw_data/sample_name
    use_vep_db - whether to use online VEP database lookups (1 for true)
          db mode a) uses online database (so cache isn't installed) b) does not use tmp files
          It is meant to be used for testing and lightweight applications.  Use the cache for
          better performance.  See discussion: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html 
    vep_cache_dir - VEP cache directory, if not doing online VEP db lookups.  Default is "/data/D_VEP"
    submit_cmd - command to initiate execution of generated script.  Default value 'bash', can set as 'cat' to allow step-by-step execution for debugging
        Replace with -d flag in future
    output_vep - write final annotated merged file in VEP rather than VCF format
    annotate_intermediate - VEP-annotate intermediate output files
    varscan_config_snp - path to varscan_snp.ini file, required for varscan run
    varscan_config_indel - path to varscan_indel.ini file, required for varscan run
    pindel_config - path to pindel.ini file, required for pindel parsing
    centromere_bed - path to BED file describing centromere regions to exclude for pindel analysis.  See SomaticWrapper/C_Centromeres
        Default: sw_dir/image.setup/C_Centromeres/pindel-centromere-exclude.bed
    gatk - path to GATK Jar file.  Default: /usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
    perl - path to PERL executable.  Default: /usr/bin/perl
    vep_cmd - path to ensembl vep executable.  Default: /usr/local/ensembl-vep/vep
    pindel_dir - path to Pindel installation dir.  Default: /usr/local/pindel
    snpsift_jar - default: /usr/local/snpEff/SnpSift.jar
    varscan_jar - default: /usr/local/VarScan.jar
    picard_jar - default: /usr/local/picard.jar
    status_rg - bam having read group or not: 1, yes and 0, no (default 1) (TODO: rename, describe in terms of pre-processing performed)
OUT

die $usage unless @ARGV >= 2;
my ( $step_number, $config_file, $config_file2 ) = @ARGV;

print("Reading configuration file $config_file\n");

# get paras from config file
# for a "key = value" pair of "xxx.yyy.zzz = foo", generates entry $params{'zzz'}='foo'
open(CONFIG, $config_file) or die "Could not open file '$config_file' $!";
my (%paras);
map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<CONFIG>);
close(CONFIG);

# Goal of an optional secondary configuration file is to allow for global and local configuration files
if (defined $config_file2) {
    if (-e $config_file2) {
        print("Reading secondary configuration file $config_file2\n");
        open(CONFIG2, $config_file2);
        map { chomp;  if(!/^[#;]/ && /=/) { @_ = split /=/; $_[1] =~ s/ //g; my $v = $_[1]; $_[0] =~ s/ //g; $paras{ (split /\./, $_[0])[-1] } = $v } } (<CONFIG2>);
        close(CONFIG2);
    } else {
        print("Optional configuration file $config_file2 does not exist.  Continuing. \n");
    }
}

map { print; print "\t"; print $paras{$_}; print "\n" } keys %paras;

# Test for mandatory parameters
die("normal_bam undefined in $config_file\n") unless exists $paras{'normal_bam'}; my $normal_bam = $paras{'normal_bam'};
die("reference_fasta undefined in $config_file\n") unless exists $paras{'reference_fasta'}; my $ref = $paras{'reference_fasta'};
die("sample_name undefined in $config_file\n") unless exists $paras{'sample_name'}; my $sample_name = $paras{'sample_name'};
# Note that some arguments which are required only for specific steps are defined when step is called, e.g. pindel_config

die("Normal BAM $normal_bam does not exist \n") if not -e $normal_bam;
die("Reference $reference_fasta does not exist \n") if not -e $reference_fasta;

# Optional arguments.  We define all default values here
my $centromere_bed="$sw_dir/image.setup/C_Centromeres/pindel-centromere-exclude.bed"; if (exists $paras{'centromere_bed'} ) { $centromere_bed=$paras{'centromere_bed'}; }
my $sw_dir="/usr/local/somaticwrapper"; if (exists $paras{'sw_dir'} ) { $sw_dir=$paras{'sw_dir'}; }
my $bsub = "bash"; if (exists $paras{'submit_cmd'} ) { $bsub=$paras{'submit_cmd'}; }
my $gatk="/usr/local/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"; if (exists $paras{'gatk'} ) { $gatk=$paras{'gatk'}; }
my $perl = "/usr/bin/perl"; if (exists $paras{'perl'} ) { $perl=$paras{'perl'}; }
my $vep_cmd="/usr/local/ensembl-vep/vep"; if (exists $paras{'vep_cmd'} ) { $vep_cmd=$paras{'vep_cmd'}; }
my $pindel_dir="/usr/local/pindel"; if (exists $paras{'pindel_dir'} ) { $pindel_dir=$paras{'pindel_dir'}; }
my $sw_data="/data/data"; if (exists $paras{'sw_data'} ) { $sw_data=$paras{'sw_data'}; }
my $snpsift_jar="/usr/local/snpEff/SnpSift.jar"; if (exists $paras{'snpsift_jar'} ) { $snpsift_jar=$paras{'snpsift_jar'}; }
my $varscan_jar="/usr/local/VarScan.jar"; if (exists $paras{'varscan_jar'} ) { $varscan_jar=$paras{'varscan_jar'}; }
my $picard="/usr/local/picard.jar"; if (exists $paras{'picard_jar'} ) { $varscan_jar=$paras{'picard_jar'}; }
my $status_rg = 1; if (exists $paras{'status_rg'} ) { $varscan_jar=$paras{'status_rg'}; }
my $gvip_dir="$sw_dir/GenomeVIP"; # GenomeVIP is not distributed separately so hard code the path

# We are retaining VEP plumbing to simplify future implementation
# die("assembly undefined in $config_file\n") unless exists $paras{'assembly'};
# my $assembly = $paras{'assembly'};
my $use_vep_db=0; if (exists $paras{'use_vep_db'} ) { $use_vep_db=$paras{'use_vep_db'}; }
my $vep_cache_dir = "/data/D_VEP"; if (exists $paras{'vep_cache_dir'} ) { $vep_cache_dir=$paras{'vep_cache_dir'}; }
my $output_vep = 0; if (exists $paras{'output_vep'} ) { $output_vep=$paras{'output_vep'}; }
my $annotate_intermediate=0; if (exists $paras{'annotate_intermediate'} ) { $annotate_intermediate=$paras{'annotate_intermediate'}; }

### begin to process each sample ###
my $sample_full_path = $sw_data."/".$sample_name;

# automatically generated scripts in runtime
my $job_files_dir="$sample_full_path/runtime";
system("mkdir -p $job_files_dir");

print("Using reference $ref\n");
print("SomaticWrapper dir: $sw_dir \n");
print("Analysis dir: $sample_full_path\n");
print("Run script dir: $job_files_dir\n");


print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
if ($step_number == 1) {
    run_GATK($normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $gatk, $picard, $status_rg);
} elsif ($step_number == 2) {
    die("varscan_config_snp undefined in $config_file\n") unless exists $paras{'varscan_config_snp'};
    die("varscan_config_indel undefined in $config_file\n") unless exists $paras{'varscan_config_indel'};
    run_varscan($normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $paras{'varscan_config_snp'}, $paras{'varscan_config_indel'});
} elsif ($step_number == 3) {
    run_pindel($normal_bam, $sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $pindel, $centromere_bed);
} elsif ($step_number == 4) {
    die("pindel_config undefined in $config_file\n") unless exists $paras{'pindel_config'};
    parse_pindel($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $perl, $gvip_dir, $paras{'pindel_config'});
} elsif ($step_number == 5) {
    merge_vcf($sample_name, $sample_full_path, $job_files_dir, $bsub, $ref, $gatk);
} else  {
    print $usage;
    die("Unknown step $step_number\n");
}


