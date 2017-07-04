## pipeline for germline_calling.pl ##
#  germline_calling.pl #

#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
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
$yellow     Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$red      	 [1]  Run gatk
$red         [2]  Run varscan
$red 		 [3]  Run Pindel
$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}
die $usage unless ($step_number >=0)&&(($step_number <= 10));
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
#my $run_script_path = `dirname $0`;
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

my $STRELKA_DIR="/gscmnt/gc2525/dinglab/rmashl/Software/bin/strelka/1.0.14/bin";
my $h37_REF="/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37V1/GRCh37-lite-chr_with_chrM.fa";
my $f_exac="/gscmnt/gc2741/ding/qgao/tools/vcf2maf-1.6.11/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz";
my $h37_REF_bai="/gscmnt/gc3027/dinglab/medseq/fasta/GRCh37/GRCh37-lite-chr_with_chrM.fa.fai";
my $pindel="/gscuser/qgao/tools/pindel/pindel";
my $PINDEL_DIR="/gscuser/qgao/tools/pindel";
#my $gatk="/gscuser/scao/tools/GenomeAnalysisTK.jar";
my $gatkexe="/gscmnt/gc2525/dinglab/rmashl/Software/bin/gatk/3.7/GenomeAnalysisTK.jar";
my $f_centromere="/gscmnt/gc3015/dinglab/medseq/Jiayin_Germline_Project/PCGP/data/pindel-centromere-exclude.bed";

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

# check to make sure the input directory has correct structure
#&check_input_dir($run_dir);
# start data processsing

if ($step_number < 5) {
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
                   #&bsub_vep();
                }
                 elsif ($step_number == 1) {
                    &bsub_gatk();
                } elsif ($step_number == 2) {
                    &bsub_varscan(1);
                }elsif ($step_number == 3) {
                    &bsub_pindel(1);
                }
				}
				}
		}
	}

sub bsub_gatk{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_gatk_g_".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    if (! -e $IN_bam_T) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_T does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_T) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    }
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
    print GATK "scr_t0=\`date \+\%s\`\n";
    print GATK "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print GATK "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
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
    print GATK "echo \"$IN_bam_N\" > \${BAMLIST}\n";
	print GATK "java  \${JAVA_OPTS} -jar "."$gatkexe -R $h37_REF"."  -T HaplotypeCaller -I \${BAMLIST} -mbq  10  -rf DuplicateRead  -rf UnmappedRead  -stand_call_conf 10.0  -o  \${rawvcf}\n";
	#print GATK "     ".$run_script_path.'genomevip_label.pl GATK \$rawvcf \$gvipvcf'."\n";
	#print GATK 'java \$JAVA_OPTS -jar '."$gatkexe -R $h37_REF".' -T SelectVariants  -V  \$gvipvcf  -o  \$snvvcf  -selectType SNP -selectType MNP'."\n";
	#print GATK 'java \$JAVA_OPTS -jar '."$gatkexe -R $h37_REF".' -T SelectVariants  -V  \$gvipvcf   -o  \$indelvcf  -selectType INDEL'."\n";
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
    if (! -e $IN_bam_T) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam file for bwa:\n";
        print "File $IN_bam_T does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    if (! -s $IN_bam_T) {#make sure input fasta file is not empty
        print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        die "Warning: Died because $IN_bam_T is empty!", $normal, "\n\n";
    }
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
    #print VARSCAN "chralt=\${chr\/:\/_}\n";
    #print VARSCAN "dir=\$chralt\n";
    print VARSCAN "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print VARSCAN "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print VARSCAN "myRUNDIR=".$sample_full_path."/varscan\n";
	print VARSCAN "outsnp=".$sample_full_path."/varscan/".$sample_name."raw.snp.vcf\n";
	print VARSCAN "logsnp=".$sample_full_path."/varscan/".$sample_name."raw.snp.log\n";
    print VARSCAN "outindel=".$sample_full_path."/varscan/".$sample_name."raw.indel.vcf\n";
    print VARSCAN "logindel=".$sample_full_path."/varscan/".$sample_name."raw.indel.log\n";
   # print VARSCAN "STATUSDIR=".$sample_full_path."/status\n";
    #print VARSCAN "RESULTSDIR=".$sample_full_path."/varscan_results\n";
    print VARSCAN "RUNDIR=".$sample_full_path."\n";
    #print VARSCAN "numgps=10\n";
    #print VARSCAN "SEQS=\"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y\"\n";
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
    #print PINDEL "echo \"$IN_bam_T\t500\t$sample_name.T\" > \${CONFIG}\n";
    print PINDEL "echo \"$IN_bam_N\t500\t$sample_name.N\" >> \${CONFIG}\n";
    print PINDEL "$pindel -T 4 -f $h37_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"." -m 6 -w 1 -J $f_centromere\n";
    close PINDEL;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
    }
 
