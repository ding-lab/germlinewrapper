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
$yellow     Usage: perl $0  --srg --step --sre --rdir --ref $normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<srg> = bam having read group or not: 1, yes and 0, no (default 1)
<sre> = re-run: 1, yes and 0, no  (default 0)
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
with chr: /gscmnt/gc3027/dinglab/medseq/fasta/GRCh37V1/GRCh37-lite-chr_with_chrM.fa
without chr: /gscmnt/gc3027/dinglab/medseq/fasta/GRCh37/GRCh37-lite.fa
mmy: /gscmnt/gc2737/ding/Reference/hs37d5_plusRibo_plusOncoViruses_plusERCC.20170530.fa 

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

my $h37_REF="";

#__PARSE COMMAND LINE
my $status = &GetOptions (
      "step=i" => \$step_number,
      "srg=i" => \$status_rg,
      "sre=i" => \$status_rerun,
      "rdir=s" => \$run_dir,
      "ref=s"  => \$h37_REF,
      "help" => \$help,
    );

#print $status,"\n";

if ($help || $run_dir eq "" || $step_number<0 || $step_number>8) {
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
my $HOME1="/gscmnt/gc2524/dinglab";
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

if($step_number==7 || $step_number==0)
    {

    print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
    $hold_job_file=$current_job_file;
    $current_job_file = "Run_report_gl_".$working_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
    print REPRUN "#!/bin/bash\n";
    print REPRUN "#BSUB -n 1\n";
    print REPRUN "#BSUB -R \"rusage[mem=40000]\"","\n";
    print REPRUN "#BSUB -M 40000000\n";
    print REPRUN "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print STREKA "#BSUB -q long\n";
    print REPRUN "#BSUB -q research-hpc\n";
    #print REPRUN "#BSUB -q ding-lab\n";
    print REPRUN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print REPRUN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print REPRUN "#BSUB -J $current_job_file\n";
    print REPRUN "#BSUB -w \"$hold_job_file\"","\n";
    print REPRUN "      ".$run_script_path."generate_final_report.pl ".$run_dir."\n";
    close REPRUN;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

}
sub bsub_gatk{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_gatk_g_".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
	my $IN_bam_N_rg=$sample_full_path."/".$sample_name.".N.rg.bam";
    #if (! -e $IN_bam_T) {#make sure there is a input fasta file 
    #    print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
    #    print "Warning: Died because there is no input bam file for bwa:\n";
    #    print "File $IN_bam_T does not exist!\n";
    #    die "Please check command line argument!", $normal, "\n\n";

    #}
   # if (! -s $IN_bam_T) {#make sure input fasta file is not empty
    #    print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
    #    die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    #}
    if (! -e $IN_bam_N) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_N does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_N) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_N is empty!", $normal, "\n\n";
    }
    open(GATK, ">$job_files_dir/$current_job_file") or die $!;
    print GATK "#!/bin/bash\n";
    print GATK "#BSUB -n 1\n";
    print GATK "#BSUB -R \"rusage[mem=30000]\"","\n";
    print GATK "#BSUB -M 30000000\n";
    print GATK "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print GATK "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print GATK "#BSUB -J $current_job_file\n";
    print GATK "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print GATK "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
  	print GATK "NBAM_rg=".$sample_full_path."/".$sample_name.".N.rg.bam\n";
	print GATK "NBAM_rg_bai=".$sample_full_path."/".$sample_name.".N.rg.bam.bai\n";
	print GATK "myRUNDIR=".$sample_full_path."/gatk\n";
    print GATK "rawvcf=".$sample_full_path."/gatk/".$sample_name.".raw.vcf\n";
	print GATK "gvipvcf=".$sample_full_path."/gatk/".$sample_name.".gvip.vcf\n";
	print GATK "snvvcf=".$sample_full_path."/gatk/".$sample_name.".snv.gvip.vcf\n";
	print GATK "indelvcf=".$sample_full_path."/gatk/".$sample_name.".indel.gvip.vcf\n";
    print GATK "RUNDIR=".$sample_full_path."\n";
    print GATK "CONFDIR="."/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/config\n";
    print GATK "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print GATK "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print GATK "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print GATK "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print GATK "if [ ! -d \${myRUNDIR} ]\n";
    print GATK "then\n";
    print GATK "mkdir \${myRUNDIR}\n";
    print GATK "fi\n";
    print GATK "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
    print GATK "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print GATK "else\n";
    print GATK "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print GATK "fi\n";
    print GATK "BAMLIST=\${RUNDIR}/gatk/bamfilelist.inp\n";
    print GATK "if [ ! -e \${BAMLIST} ]\n";
    print GATK "then\n";
    print GATK "rm \${BAMLIST}\n";
    print GATK "fi\n";
   # print GATK "echo \"$IN_bam_N_rg\" > \${BAMLIST}\n";
    print GATK "if [ $status_rg -eq 0 ]\n";
	print GATK "then\n";
	print GATK "java  \${JAVA_OPTS} -jar "."$picardexe AddOrReplaceReadGroups I=\${NBAM} O=\${NBAM_rg} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20\n";
	print GATK "samtools index \${NBAM_rg}\n";
   # print GATK "else\n";
	print GATK "java  \${JAVA_OPTS} -jar "."$gatkexe -R $h37_REF"."  -T HaplotypeCaller -I \${NBAM_rg} -mbq  10  -rf DuplicateRead  -rf UnmappedRead  -stand_call_conf 10.0  -o  \${rawvcf}\n";
    print GATK "rm \${NBAM_rg}\n";
    print GATK "rm \${NBAM_rg_bai}\n";
	print GATK "else\n";
	print GATK "java  \${JAVA_OPTS} -jar "."$gatkexe -R $h37_REF"."  -T HaplotypeCaller -I \${NBAM} -mbq  10  -rf DuplicateRead  -rf UnmappedRead  -stand_call_conf 10.0  -o  \${rawvcf}\n";
    print GATK "fi\n";
	print GATK "     ".$run_script_path."genomevip_label.pl GATK \${rawvcf} \${gvipvcf}"."\n";
	print GATK "java \${JAVA_OPTS} -jar "."$gatkexe -R $h37_REF"." -T SelectVariants  -V  \${gvipvcf}  -o  \${snvvcf}  -selectType SNP -selectType MNP"."\n";
	print GATK "java \${JAVA_OPTS} -jar "."$gatkexe -R $h37_REF"." -T SelectVariants  -V  \${gvipvcf}   -o  \${indelvcf}  -selectType INDEL"."\n";	
	close GATK;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

	}

sub bsub_varscan{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j2_varscan_g_".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    #if (! -e $IN_bam_T) {#make sure there is a input fasta file 
    #    print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
    #    print "Warning: Died because there is no input bam file for bwa:\n";
    #    print "File $IN_bam_T does not exist!\n";
    #    die "Please check command line argument!", $normal, "\n\n";

   # }
   # if (! -s $IN_bam_T) {#make sure input fasta file is not empty
    #    print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
     #   die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    #}
    if (! -e $IN_bam_N) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_N does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_N) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_N is empty!", $normal, "\n\n";
    }
    open(VARSCAN, ">$job_files_dir/$current_job_file") or die $!;
    print VARSCAN "#!/bin/bash\n";
    print VARSCAN "#BSUB -n 1\n";
    print VARSCAN "#BSUB -R \"rusage[mem=30000]\"","\n";
    print VARSCAN "#BSUB -M 30000000\n";
    print VARSCAN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print VARSCAN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print VARSCAN "#BSUB -J $current_job_file\n";
    print VARSCAN "#BSUB -w \"$hold_job_file\"","\n";
    print VARSCAN "scr_t0=\`date \+\%s\`\n";
    print VARSCAN "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print VARSCAN "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print VARSCAN "myRUNDIR=".$sample_full_path."/varscan\n";
	print VARSCAN "outsnp=".$sample_full_path."/varscan/".$sample_name."raw.snp.vcf\n";
	print VARSCAN "logsnp=".$sample_full_path."/varscan/".$sample_name."raw.snp.log\n";
    print VARSCAN "outindel=".$sample_full_path."/varscan/".$sample_name."raw.indel.vcf\n";
    print VARSCAN "logindel=".$sample_full_path."/varscan/".$sample_name."raw.indel.log\n";
    print VARSCAN "RUNDIR=".$sample_full_path."\n";
    print VARSCAN "CONFDIR="."/gscmnt/gc2521/dinglab/cptac_prospective_samples/exome/config\n";
    print VARSCAN "GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip\n";
    print VARSCAN "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print VARSCAN "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print VARSCAN "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_60-x64\n";
    print VARSCAN "export JAVA_OPTS=\"-Xms256m -Xmx512m\"\n";
    print VARSCAN "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print VARSCAN "if [ ! -d \${myRUNDIR} ]\n";
    print VARSCAN "then\n";
    print VARSCAN "mkdir \${myRUNDIR}\n";
    print VARSCAN "fi\n";
    print VARSCAN "if \[\[ -z \"\$LD_LIBRARY_PATH\" \]\] \; then\n";
    print VARSCAN "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib\n";
    print VARSCAN "else\n";
    print VARSCAN "export LD_LIBRARY_PATH=\${JAVA_HOME}/lib:\${LD_LIBRARY_PATH}\n";
    print VARSCAN "fi\n";
    print VARSCAN "BAMLIST=\${RUNDIR}/varscan/bamfilelist.inp\n";
    print VARSCAN "if [ ! -e \${BAMLIST} ]\n";
    print VARSCAN "then\n";
    print VARSCAN "rm \${BAMLIST}\n";
    print VARSCAN "fi\n";
    print VARSCAN "echo \"$IN_bam_N\" > \${BAMLIST}\n";
	print VARSCAN "ncols=\$(echo \"3*( \$(wc -l < \$BAMLIST) +1)\"|bc)\n";
    print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar mpileup2snp  -  --p-value  0.10   --min-coverage  3   --min-var-freq  0.08   --min-reads2  2   --min-avg-qual  15   --min-freq-for-hom  0.75   --strand-filter  1   --output-vcf  1   > \${outsnp}  2> \${logsnp}\n";   
 	print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar mpileup2indel  -  --p-value  0.10   --min-coverage  3   --min-var-freq  0.20   --min-reads2  2   --min-avg-qual  15   --min-freq-for-hom  0.75   --strand-filter  1   --output-vcf  1   > \${outindel}  2> \${logindel}\n";
    close VARSCAN;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

	}

sub bsub_pindel{
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j3_pindel_g_".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    open(PINDEL, ">$job_files_dir/$current_job_file") or die $!;
    print PINDEL "#!/bin/bash\n";
    print PINDEL "#BSUB -n 4\n";
    print PINDEL "#BSUB -R \"span[hosts=1] rusage[mem=30000]\"","\n";
    print PINDEL "#BSUB -M 30000000\n";
    print PINDEL "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print PINDEL "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print PINDEL "#BSUB -J $current_job_file\n";
    print PINDEL "#BSUB -w \"$hold_job_file\"","\n";
    print PINDEL "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print PINDEL "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print PINDEL "myRUNDIR=".$sample_full_path."/pindel\n";
    print PINDEL "CONFIG=\${myRUNDIR}"."/".$sample_name.".config\n";
    print PINDEL "if [ ! -d \${myRUNDIR} ]\n";
    print PINDEL "then\n";
    print PINDEL "mkdir \${myRUNDIR}\n";
    print PINDEL "fi\n";
    print PINDEL "echo \"$IN_bam_N\t500\t$sample_name.N\" >> \${CONFIG}\n";
    print PINDEL "$pindel -T 4 -f $h37_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"." -m 6 -w 1 -J $f_centromere\n";
    close PINDEL;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
    }

sub bsub_parse_pindel {

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j4_parse_pindel_g_".$sample_name.".sh";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    open(PP, ">$job_files_dir/$current_job_file") or die $!;
    print PP "#!/bin/bash\n";
    print PP "#BSUB -n 1\n";
    print PP "#BSUB -R \"rusage[mem=30000]\"","\n";
    print PP "#BSUB -M 30000000\n";
    print PP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print PP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print PP "#BSUB -J $current_job_file\n";
    #print PP "#BSUB -q long\n";
    print PP "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    print PP "#BSUB -q research-hpc\n";
    print PP "#BSUB -w \"$hold_job_file\"","\n";
    print PP "RUNDIR=".$sample_full_path."\n";
    print PP "cat > \${RUNDIR}/pindel/pindel_filter.input <<EOF\n";
    print PP "pindel.filter.pindel2vcf = $PINDEL_DIR/pindel2vcf\n";
    print PP "pindel.filter.variants_file = \${RUNDIR}/pindel/pindel.out.raw\n";
    print PP "pindel.filter.REF = $h37_REF\n";
    print PP "pindel.filter.date = 000000\n";
    print PP "pindel.filter.heterozyg_min_var_allele_freq = 0.2\n";
    print PP "pindel.filter.homozyg_min_var_allele_freq = 0.8\n";
    print PP "pindel.filter.mode = germline\n";
    print PP "pindel.filter.apply_filter = true\n";
    print PP "pindel.filter.germline.min_coverages = 10\n";
    print PP "pindel.filter.germline.min_var_allele_freq = 0.20\n";
    print PP "pindel.filter.germline.require_balanced_reads = \"true\"\n";
    print PP "pindel.filter.germline.remove_complex_indels = \"true\"\n";
    print PP "pindel.filter.germline.max_num_homopolymer_repeat_units = 6\n";
    print PP "EOF\n";
    print PP "myRUNDIR=".$sample_full_path."/pindel\n";
    print PP "cd \${RUNDIR}/pindel\n";
    print PP "outlist=pindel.out.filelist\n";
    print PP "find \. -name \'*_D\' -o -name \'*_SI\' -o -name \'*_INV\' -o -name \'*_TD\'  > \./\${outlist}\n";
    print PP 'list=$(xargs -a  ./$outlist)'."\n";
    print PP "pin_var_file=pindel.out.raw\n";
    print PP 'cat $list | grep ChrID > ./$pin_var_file'."\n";
    print PP "     ".$run_script_path."pindel_filter.v0.5.pl ./pindel_filter.input\n";
    close PP;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";		
    system ($bsub_com);

    }


sub bsub_merge_vcf{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j5_merge_vcf_g.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    print MERGE "#BSUB -n 1\n";
    print MERGE "#BSUB -R \"rusage[mem=30000]\"","\n";
    print MERGE "#BSUB -M 30000000\n";
    print MERGE "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print MERGE "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print MERGE "#BSUB -J $current_job_file\n";
    print MERGE "#BSUB -q long\n";
   # print MERGE "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    #print MERGE "#BSUB -q research-hpc\n";
    print MERGE "#BSUB -w \"$hold_job_file\"","\n";
  	print MERGE "RUNDIR=".$sample_full_path."\n";
    #print VEP "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print MERGE "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print MERGE "export JAVA_HOME=/gscmnt/gc2525/dinglab/rmashl/Software/bin/jre/1.8.0_121-x64\n";
    print MERGE "export JAVA_OPTS=\"-Xmx10g\"\n";
    print MERGE "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print MERGE "GATK_snv_VCF="."\${RUNDIR}/gatk/$sample_name.snv.gvip.filtered.vcf\n";
    print MERGE "GATK_indel_VCF="."\${RUNDIR}/gatk/$sample_name.indel.gvip.filtered.vcf\n";
    print MERGE "VARSCAN_snv_VCF="."\${RUNDIR}/varscan/".$sample_name."raw.snp.filtered.vcf\n";
	print MERGE "VARSCAN_indel_VCF="."\${RUNDIR}/varscan/".$sample_name."raw.indel.filtered.vcf\n";	
    print MERGE "PINDEL_VCF="."\${RUNDIR}/pindel/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf\n";
    print MERGE "MERGER_OUT="."\${RUNDIR}/merged.vcf\n";
    print MERGE "cat > \${RUNDIR}/vep.merged.input <<EOF\n";
    print MERGE "merged.vep.vcf = ./merged.filtered.vcf\n";
    print MERGE "merged.vep.output = ./merged.VEP.vcf\n";
    print MERGE "merged.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print MERGE "merged.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print MERGE "merged.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print MERGE "merged.vep.assembly = GRCh37\n";
    print MERGE "EOF\n";
   # print MERGE "java \${JAVA_OPTS} -jar $gatk -R $h37_REF -T CombineVariants -o \${MERGER_OUT} --variant:gsnp \${GATK_snv_VCF} --variant:gindel \${GATK_indel_VCF} --variant:vsnp \${VARSCAN_snv_VCF} --variant:vindel \${VARSCAN_indel_VCF} --variant:pindel \${PINDEL_VCF} -genotypeMergeOptions UNIQUIFY\n"; 
    print MERGE "     ".$run_script_path."filter_gatk_varscan.pl \${RUNDIR} $sample_name\n";
	print MERGE "java \${JAVA_OPTS} -jar $gatk -R $h37_REF -T CombineVariants -o \${MERGER_OUT} --variant:gsnp \${GATK_snv_VCF} --variant:gindel \${GATK_indel_VCF} --variant:vsnp \${VARSCAN_snv_VCF} --variant:vindel \${VARSCAN_indel_VCF} --variant:pindel \${PINDEL_VCF} -genotypeMergeOptions PRIORITIZE -priority gsnp,vsnp,gindel,vindel,pindel\n";	
#-priority gsnp,vsnp,gindel,vindel,pindel\n";
    #print MERGE "     ".$run_script_path."vaf_filter.pl \${RUNDIR}\n";
    #print MERGE "cd \${RUNDIR}\n";
    #print MERGE ". /gscmnt/gc2525/dinglab/rmashl/Software/perl/set_envvars\n";
    #print MERGE "     ".$run_script_path."vep_annotator.pl ./vep.merged.input >&./vep.merged.log\n";
    close MERGE;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #$bsub_com = "sh $job_files_dir/$current_job_file\n";
    system ($bsub_com);
    
	}


sub bsub_vcf_2_maf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j6_vcf_2_maf.".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    #my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    open(MAF, ">$job_files_dir/$current_job_file") or die $!;
    print MAF "#!/bin/bash\n";
    print MAF "#BSUB -n 1\n";
    print MAF "#BSUB -R \"rusage[mem=30000]\"","\n";
    print MAF "#BSUB -M 30000000\n";
    print MAF "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print MAF "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print MAF "#BSUB -J $current_job_file\n";
    print MAF "#BSUB -q long\n";
    print MAF "#BSUB -w \"$hold_job_file\"","\n";
    print MAF "RUNDIR=".$sample_full_path."\n";
    print MAF "cat > \${RUNDIR}/vep.merged.input <<EOF\n";
    print MAF "merged.vep.vcf = ./merged.vcf\n";
    print MAF "merged.vep.output = ./merged.VEP.vcf\n";
    print MAF "merged.vep.vep_cmd = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print MAF "merged.vep.cachedir = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache\n";
    print MAF "merged.vep.reffasta = /gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa\n";
    print MAF "merged.vep.assembly = GRCh37\n";
    print MAF "EOF\n";
    print MAF "F_VCF_1=".$sample_full_path."/merged.vcf\n";
    print MAF "F_VCF_2=".$sample_full_path."/".$sample_name.".vcf\n";
    print MAF "F_VEP_1=".$sample_full_path."/merged.VEP.vcf\n";
    print MAF "F_VEP_2=".$sample_full_path."/".$sample_name.".vep.vcf\n";
    print MAF "F_maf=".$sample_full_path."/".$sample_name.".maf\n";
    print MAF "cd \${RUNDIR}\n";
    print MAF ". /gscmnt/gc2525/dinglab/rmashl/Software/perl/set_envvars\n";
    print MAF "     ".$run_script_path."vep_annotator.pl ./vep.merged.input >&./vep.merged.log\n";
    print MAF "rm \${F_VCF_2}\n";
    print MAF "rm \${F_VEP_2}\n";
    print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
    print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
    print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $h37_REF --filter-vcf $f_exac\n";
    close MAF;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

    }
