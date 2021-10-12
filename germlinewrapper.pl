#  the input bam without readgroup information #	
### Song Cao ###
### last updated 7/6/2017 ###
### last updated 10/09/2017 ###
### last updated 10/19/2017 ###
### updated 1/9/2018 ###
## cp1, 10/2011 ##
#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;
my $version = 1.1;
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
$yellow     Usage: perl $0  --srg --step --sre --rdir --ref --log --q

$normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<srg> = bam having read group or not: 1, yes and 0, no (default 1)
<log> = full path of the folder for saving log file; usually upper folder of rdir 
<sre> = re-run: 1, yes and 0, no  (default 0)
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
<q> which queue for submitting job; research-hpc, ding-lab, long (default)

GDC HG38: /storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa 

<run_folder> = full path of the folder holding files for this sequence run
<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$red      	 [1]  Run gatk
$red         	 [2]  Run varscan
$red 		 [3]  Run pindel
$yellow 	 [4]  Parse pindel
$yellow          [5]  filter vcf
$purple 	 [6]  Merge calls
$green 		 [7]  VCF2MAF
$cyan 		 [8]  Generate final maf

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
my $h38_REF="";
my $q_name="";
my $chr_status=0; 

my $compute_username="songcao";
my $group_name="testing";

my $status = &GetOptions (
      "step=i" => \$step_number,
      "srg=i" => \$status_rg,
      "sre=i" => \$status_rerun,
      "rdir=s" => \$run_dir,
      "ref=s"  => \$h38_REF,
      "log=s"  => \$log_dir,
	  "q=s" => \$q_name,
      "help" => \$help,
    );

if ($help || $run_dir eq "" || $log_dir eq ""  || $step_number<0 || $step_number>8) {
      print $usage;
      exit;
   }

print "run dir=",$run_dir,"\n";
print "step num=",$step_number,"\n";
print "status rerun=",$status_rerun,"\n";
print "status readgroup=",$status_rg,"\n";
print "queue name=",$q_name,"\n";
if($q_name eq "")
{
    $q_name="long";
}

if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

my $email = "scao\@wustl\.edu";
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-2];
my $HOME1=$log_dir;

if (! -d $HOME1)
{
`mkdir $HOME1`;
}

if (! -d $HOME1."/tmpgermline") {
    `mkdir $HOME1"/tmpgermline"`;
}

my $job_files_dir = $HOME1."/tmpgermline";

if (! -d $HOME1."/LSF_DIR_GERMLINE") {
    `mkdir $HOME1"/LSF_DIR_GERMLINE"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR_GERMLINE";

my $run_script_path =`echo \$PWD`;
chomp $run_script_path;
my $script_dir=$run_script_path;

$run_script_path = "/usr/bin/perl ".$run_script_path."/";
print $run_script_path,"\n";
my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
my $h38_REF_bai=$h38_REF.".fai";
my $gatk="/storage1/fs1/songcao/Active/Software/GenomeAnalysis/GenomeAnalysisTK.jar";
my $samtools="/storage1/fs1/songcao/Active/Software/anaconda3/bin";
my $varscan="/storage1/fs1/songcao/Active/Software/varscan/2.3.8";
my $STRELKA_DIR="/gscmnt/gc2525/dinglab/rmashl/Software/bin/strelka/1.0.14/bin";
my $f_ref_annot="/gscmnt/gc2518/dinglab/scao/tools/vep/Homo_sapiens.GRCh38.dna.primary_assembly.fa";
my $vepcache="/gscmnt/gc2518/dinglab/scao/tools/vep/v85";
my $pindel="/storage1/fs1/songcao/Active/Software/anaconda3/bin/pindel";
my $PINDEL_DIR="/storage1/fs1/songcao/Active/Software/anaconda3/bin";
my $gatkexe3="/storage1/fs1/songcao/Active/Software/gatk/3.7/GenomeAnalysisTK.jar";
my $gatkexe4="gatk";
my $picardexe="/storage1/fs1/songcao/Active/Software/picard/picard.jar";
my $java_dir="/storage1/fs1/songcao/Active/Software/jre1.8.0_121";

my $vepcmd="/gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v85/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl";

my $first_line=`head -n 1 $h38_REF`; 

if($first_line=~/^\>chr/) { $chr_status=1; }

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

if ($step_number < 8) {
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
				   &bsub_filter_vcf();
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
                    &bsub_filter_vcf(1);
                }elsif ($step_number == 6) {
                    &bsub_merge_vcf(1);
                }elsif ($step_number==7) {
				    &bsub_vcf_2_maf(1);		
				}
				}
				}
		}
	}

if($step_number==8 || $step_number==0)
    {

    print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
    $hold_job_file=$current_job_file;
    $current_job_file = "Run_report_gl_".$working_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    if(-e $lsf_out) 
	{
	`rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
	}

    open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
    print REPRUN "#!/bin/bash\n";
    #print REPRUN "#BSUB -n 1\n";
    #print REPRUN "#BSUB -R \"rusage[mem=40000]\"","\n";
    #print REPRUN "#BSUB -M 40000000\n";
    #print REPRUN "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print REPRUN "#BSUB -q long\n";
    #print REPRUN "#BSUB -q research-hpc\n";
    #print REPRUN "#BSUB -q ding-lab\n";
    #print REPRUN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print REPRUN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print REPRUN "#BSUB -J $current_job_file\n";
    #print REPRUN "#BSUB -w \"$hold_job_file\"","\n";
    print REPRUN "      ".$run_script_path."generate_final_report.pl ".$run_dir."\n";
    close REPRUN;
    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ($bsub_com);

 	my $sh_file=$job_files_dir."/".$current_job_file;

    if($q_name eq "research-hpc")
    {
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";     }
    else {        $bsub_com = "bsub -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";   }
    print $bsub_com;
    system ($bsub_com);

}
sub bsub_gatk{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";

	my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
  
    $current_job_file = "j1_gatk_g_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }

    #my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
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
    #print GATK "#BSUB -n 1\n";
    #print GATK "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print GATK "#BSUB -M 30000000\n";
    #print GATK "#BSUB -q ding-lab\n";
	#print GATK "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print GATK "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print GATK "#BSUB -J $current_job_file\n";
    print GATK "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print GATK "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print GATK "NBAM_rg=".$sample_full_path."/".$sample_name.".N.rg.bam\n";
    print GATK "NBAM_rg_bai=".$sample_full_path."/".$sample_name.".N.rg.bam.bai\n";
    print GATK "myRUNDIR=".$sample_full_path."/gatk\n";
    print GATK "RUNDIR=".$sample_full_path."\n";
    print GATK "export SAMTOOLS_DIR=$samtools\n";
    print GATK "export JAVA_HOME=$java_dir\n";
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
    print GATK "if [ $status_rg -eq 0 ]\n";
    print GATK "then\n";
    print GATK "java  \${JAVA_OPTS} -jar "."$picardexe AddOrReplaceReadGroups I=\${NBAM} O=\${NBAM_rg} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20\n";
    print GATK "samtools index \${NBAM_rg}\n";
    foreach my $chr (@chrlist)
    {
	my $chr1=$chr; 
	if($chr_status==1) { $chr1="chr".$chr; }
    	print GATK "rawvcf=".$sample_full_path."/gatk/".$sample_name.".raw.$chr.vcf\n";
	print GATK "$gatkexe4 HaplotypeCaller  -I \${NBAM_rg} -L $chr1 -O \${rawvcf} -R $h38_REF -RF NotDuplicateReadFilter -RF MappingQualityReadFilter -RF MappedReadFilter\n";
	}
        print GATK "rm \${NBAM_rg}\n";
        print GATK "rm \${NBAM_rg_bai}\n";
	print GATK "else\n";
	print GATK "echo \"run gatk4\"","\n";
	foreach my $chr (@chrlist)
	{
	 my $chr1=$chr;
         if($chr_status==1) { $chr1="chr".$chr; }
         print GATK "rawvcf=".$sample_full_path."/gatk/".$sample_name.".raw.$chr.vcf\n";
	 print GATK "$gatkexe4 HaplotypeCaller  -I \${NBAM} -O \${rawvcf} -R $h38_REF -L $chr1 -RF NotDuplicateReadFilter -RF MappingQualityReadFilter -RF MappedReadFilter\n";
	}
    	print GATK "fi\n";

    	foreach my $chr (@chrlist)
    	{
    	print GATK "rawvcf=".$sample_full_path."/gatk/".$sample_name.".raw.$chr.vcf\n";
    	print GATK "gvipvcf=".$sample_full_path."/gatk/".$sample_name.".gvip.$chr.vcf\n";
    	print GATK "snvvcf=".$sample_full_path."/gatk/".$sample_name.".snv.gvip.$chr.vcf\n";
    	print GATK "indelvcf=".$sample_full_path."/gatk/".$sample_name.".indel.gvip.$chr.vcf\n";
	print GATK "     ".$run_script_path."genomevip_label.pl GATK \${rawvcf} \${gvipvcf}"."\n";	
 	print GATK "$gatkexe4 SelectVariants -R $h38_REF -V  \${gvipvcf}  -O  \${snvvcf}  -select-type SNP -select-type MNP"."\n";    
	print GATK "$gatkexe4 SelectVariants -R $h38_REF -V  \${gvipvcf}  -O  \${indelvcf}  -select-type INDEL"."\n";
	}

	print GATK "     ".$run_script_path."merge_gatk.pl $sample_full_path $sample_name\n"; 
	close GATK;

 	my $sh_file=$job_files_dir."/".$current_job_file;
	

    $bsub_com = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(broadinstitute/gatk)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

	}

sub bsub_varscan{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j2_varscan_g_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }
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
    #print VARSCAN "#BSUB -n 1\n";
    #print VARSCAN "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print VARSCAN "#BSUB -M 30000000\n";
    #print VARSCAN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print VARSCAN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
   # print VARSCAN "#BSUB -J $current_job_file\n";
    #print VARSCAN "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print VARSCAN "#BSUB -q research-hpc\n";
    #print VARSCAN "#BSUB -w \"$hold_job_file\"","\n";
    print VARSCAN "scr_t0=\`date \+\%s\`\n";
    print VARSCAN "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print VARSCAN "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print VARSCAN "myRUNDIR=".$sample_full_path."/varscan\n";
    print VARSCAN "outsnp=".$sample_full_path."/varscan/".$sample_name."raw.snp.vcf\n";
    print VARSCAN "logsnp=".$sample_full_path."/varscan/".$sample_name."raw.snp.log\n";
    print VARSCAN "outindel=".$sample_full_path."/varscan/".$sample_name."raw.indel.vcf\n";
    print VARSCAN "logindel=".$sample_full_path."/varscan/".$sample_name."raw.indel.log\n";
    print VARSCAN "RUNDIR=".$sample_full_path."\n";
    print VARSCAN "export VARSCAN_DIR=$varscan\n";
    print VARSCAN "export SAMTOOLS_DIR=$samtools\n";
    print VARSCAN "export JAVA_HOME=$java_dir\n";
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
    print VARSCAN "ncols=6\n";
    print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h38_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar mpileup2snp  -  --p-value  0.10   --min-coverage  3   --min-var-freq  0.08   --min-reads2  2   --min-avg-qual  15   --min-freq-for-hom  0.75   --strand-filter  1   --output-vcf  1   > \${outsnp}  2> \${logsnp}\n";   
    print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h38_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar mpileup2indel  -  --p-value  0.10   --min-coverage  3   --min-var-freq  0.20   --min-reads2  2   --min-avg-qual  15   --min-freq-for-hom  0.75   --strand-filter  1   --output-vcf  1   > \${outindel}  2> \${logindel}\n";
    close VARSCAN;

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);

}

sub bsub_pindel{
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j3_pindel_g_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    open(PINDEL, ">$job_files_dir/$current_job_file") or die $!;
    print PINDEL "#!/bin/bash\n";
    print PINDEL "TBAM=".$sample_full_path."/".$sample_name.".T.bam\n";
    print PINDEL "NBAM=".$sample_full_path."/".$sample_name.".N.bam\n";
    print PINDEL "myRUNDIR=".$sample_full_path."/pindel\n";
    print PINDEL "CONFIG=\${myRUNDIR}"."/".$sample_name.".config\n";
    print PINDEL "if [ ! -d \${myRUNDIR} ]\n";
    print PINDEL "then\n";
    print PINDEL "mkdir \${myRUNDIR}\n";
    print PINDEL "fi\n";
    print PINDEL "echo \"$IN_bam_N\t500\t$sample_name.N\" > \${CONFIG}\n";
    print PINDEL "$pindel -T 4 -f $h38_REF -i \${CONFIG} -o \${myRUNDIR}"."/$sample_name"." -m 6 -w 1\n";
    close PINDEL;
    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);
    
    }

sub bsub_parse_pindel {

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }

    $current_job_file = "j4_parse_pindel_g_".$sample_name.".sh";

    #my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    #my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    #`rm $lsf_out`;
    #`rm $lsf_err`;
    open(PP, ">$job_files_dir/$current_job_file") or die $!;
    print PP "#!/bin/bash\n";
    #print PP "#BSUB -n 1\n";
    #print PP "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print PP "#BSUB -M 30000000\n";
    #print PP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print PP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print PP "#BSUB -J $current_job_file\n";
    #print PP "#BSUB -q ding-lab\n";
    #print PP "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
#print PP "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/lucid-default:latest)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    #print PP "#BSUB -q research-hpc\n";
    print PP "#BSUB -w \"$hold_job_file\"","\n";
    print PP "RUNDIR=".$sample_full_path."\n";
    print PP "cat > \${RUNDIR}/pindel/pindel_filter.input <<EOF\n";
    print PP "pindel.filter.pindel2vcf = $PINDEL_DIR/pindel2vcf\n";
    print PP "pindel.filter.variants_file = \${RUNDIR}/pindel/pindel.out.raw\n";
    print PP "pindel.filter.REF = $h38_REF\n";
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

    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);
    }


sub bsub_filter_vcf{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j5_filter_vcf_g.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";


    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

   if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }

    open(FILTER, ">$job_files_dir/$current_job_file") or die $!;
    print FILTER "#!/bin/bash\n";
    #print MERGE "#BSUB -n 1\n";
    ##print MERGE "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print MERGE "#BSUB -M 30000000\n";
    #print MERGE "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print MERGE "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print MERGE "#BSUB -J $current_job_file\n";
    #print MERGE "#BSUB -q long\n";
    #print MERGE "#BSUB -q ding-lab\n"; 
    #print MERGE "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    #print MERGE "#BSUB -q research-hpc\n";
    print FILTER "#BSUB -w \"$hold_job_file\"","\n";
    print FILTER "RUNDIR=".$sample_full_path."\n";
    #print VEP "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print FILTER "export SAMTOOLS_DIR=$samtools\n";
    print FILTER "export JAVA_HOME=$java_dir\n";
    print FILTER "export JAVA_OPTS=\"-Xmx10g\"\n";
    print FILTER "export PATH=\${JAVA_HOME}/bin:\${PATH}\n";
    print FILTER "GATK_snv_VCF="."\${RUNDIR}/gatk/$sample_name.snv.gvip.filtered.vcf\n";
    print FILTER "GATK_indel_VCF="."\${RUNDIR}/gatk/$sample_name.indel.gvip.filtered.vcf\n";
    print FILTER "VARSCAN_snv_VCF="."\${RUNDIR}/varscan/".$sample_name."raw.snp.filtered.vcf\n";
    print FILTER "VARSCAN_indel_VCF="."\${RUNDIR}/varscan/".$sample_name."raw.indel.filtered.vcf\n";
    print FILTER "PINDEL_VCF="."\${RUNDIR}/pindel/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf\n";
    print FILTER "     ".$run_script_path."filter_gatk_varscan.pl \${RUNDIR} $sample_name\n";
    close FILTER;


    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>80000] rusage[mem=80000]\" -M 80000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
    print $bsub_com;
    system ($bsub_com);
}

sub bsub_merge_vcf{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j6_merge_vcf_g.".$sample_name.".sh";
    my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";


    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    
   if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }

    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    #print MERGE "#BSUB -n 1\n";
    ##print MERGE "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print MERGE "#BSUB -M 30000000\n";
    #print MERGE "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print MERGE "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print MERGE "#BSUB -J $current_job_file\n";
    #print MERGE "#BSUB -q long\n";
    #print MERGE "#BSUB -q ding-lab\n"; 
    #print MERGE "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    #print MERGE "#BSUB -q research-hpc\n";
    print MERGE "#BSUB -w \"$hold_job_file\"","\n";
  	print MERGE "RUNDIR=".$sample_full_path."\n";
    #print VEP "export VARSCAN_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/varscan/2.3.8\n";
    print MERGE "export SAMTOOLS_DIR=/gscmnt/gc2525/dinglab/rmashl/Software/bin/samtools/1.2/bin\n";
    print MERGE "export JAVA_HOME=$java_dir\n";
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
    print MERGE "merged.vep.vep_cmd = $vepcmd\n";
#="/gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v85/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_pre/gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v81/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print MERGE "merged.vep.cachedir = $vepcache\n";
    print MERGE "merged.vep.reffasta = $f_ref_annot\n";
    print MERGE "merged.vep.assembly = GRCh38\n";
    print MERGE "EOF\n";
   # print MERGE "java \${JAVA_OPTS} -jar $gatk -R $h38_REF -T CombineVariants -o \${MERGER_OUT} --variant:gsnp \${GATK_snv_VCF} --variant:gindel \${GATK_indel_VCF} --variant:vsnp \${VARSCAN_snv_VCF} --variant:vindel \${VARSCAN_indel_VCF} --variant:pindel \${PINDEL_VCF} -genotypeMergeOptions UNIQUIFY\n"; 
	### snv 
    #print MERGE "     ".$run_script_path."filter_gatk_varscan.pl \${RUNDIR} $sample_name\n";
	print MERGE "java \${JAVA_OPTS} -jar $gatk -R $h38_REF -T CombineVariants -o \${MERGER_OUT} --variant:gsnp \${GATK_snv_VCF} --variant:gindel \${GATK_indel_VCF} --variant:vsnp \${VARSCAN_snv_VCF} --variant:vindel \${VARSCAN_indel_VCF} --variant:pindel \${PINDEL_VCF} -genotypeMergeOptions PRIORITIZE -priority gsnp,vsnp,gindel,vindel,pindel\n";	
#-priority gsnp,vsnp,gindel,vindel,pindel\n";
    #print MERGE "     ".$run_script_path."vaf_filter.pl \${RUNDIR}\n";
    #print MERGE "cd \${RUNDIR}\n";
    #print MERGE ". /gscmnt/gc2525/dinglab/rmashl/Software/perl/set_envvars\n";
    #print MERGE "     ".$run_script_path."vep_annotator.pl ./vep.merged.input >&./vep.merged.log\n";
    close MERGE;
    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #$bsub_com = "sh $job_files_dir/$current_job_file\n";
    #system ($bsub_com);
 
    my $sh_file=$job_files_dir."/".$current_job_file;

    if($q_name eq "research-hpc")
    {
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>80000] rusage[mem=80000]\" -M 80000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";     }
    else {        $bsub_com = "bsub -q $q_name -n 1 -R \"select[mem>80000] rusage[mem=80000]\" -M 80000000 -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";
    }

    print $bsub_com;

    system ($bsub_com);

	}


sub bsub_vcf_2_maf{
  
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }


    $current_job_file = "j7_vcf_2_maf.".$sample_name.".sh";
    #my $IN_bam_T = $sample_full_path."/".$sample_name.".T.bam";
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out)
    {
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    }

    open(MAF, ">$job_files_dir/$current_job_file") or die $!;
    print MAF "#!/bin/bash\n";
    #print MAF "#BSUB -n 1\n";
    #print MAF "#BSUB -R \"rusage[mem=60000]\"","\n";
    #print MAF "#BSUB -M 60000000\n";
    #print MAF "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print VARSCANP "#BSUB -q long\n";
    #print MAF "#BSUB -q research-hpc\n";
    #print MAF "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print MAF "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print MAF "#BSUB -J $current_job_file\n";
    #print MAF "#BSUB -q ding-lab\n";
    #print MAF "#BSUB -w \"$hold_job_file\"","\n";
    print MAF "RUNDIR=".$sample_full_path."\n";
    print MAF "cat > \${RUNDIR}/vep.merged.input <<EOF\n";
    print MAF "merged.vep.vcf = ./merged.vcf\n";
    print MAF "merged.vep.output = ./merged.VEP.vcf\n";
    print MAF "merged.vep.vep_cmd = $vepcmd\n";
#/gscmnt/gc2525/dinglab/rmashl/Software/bin/VEP/v85/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl\n";
    print MAF "merged.vep.cachedir = $vepcache\n";
   	print MAF "merged.vep.reffasta = $f_ref_annot\n";
    print MAF "merged.vep.assembly = GRCh38\n";
    print MAF "EOF\n";
	print MAF "F_VCF_0=".$sample_full_path."/merged.vcf\n";
    print MAF "F_VCF_1=".$sample_full_path."/merged.1.vcf\n";
    print MAF "F_VCF_2=".$sample_full_path."/".$sample_name.".vcf\n";
    print MAF "F_VEP_1=".$sample_full_path."/merged.VEP.vcf\n";
    print MAF "F_VEP_2=".$sample_full_path."/".$sample_name.".vep.vcf\n";
    print MAF "F_maf=".$sample_full_path."/".$sample_name.".maf\n";
	print MAF "vep_log=".$sample_full_path."/vep.merged.log\n";
    print MAF "if [ $status_rerun -eq 1 ]\n";
    print MAF "then\n";
    print MAF "rm \${vep_log}\n";
    print MAF "fi\n";
	
    print MAF "if [ -f \${vep_log} ]\n";
	print MAF "then\n";
    print MAF 'tail -1 ${vep_log} | grep ERROR',"\n";
	print MAF '          CHECK=$?',"\n";
	print MAF '		if [ ${CHECK} -eq 0 ]',"\n";
	print MAF "then\n"; 
    print MAF "     ".$run_script_path."remove_svtype.pl \${F_VCF_0} \${F_VCF_1}\n";
    print MAF "cd \${RUNDIR}\n";
	print MAF ". $script_dir/set_envvars\n";
    print MAF "     ".$run_script_path."vep_annotator_v1.1.pl ./vep.merged.input >&./vep.merged.log\n";
    print MAF "rm \${F_VCF_2}\n";
    print MAF "rm \${F_VEP_2}\n";
    print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
    print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
#    print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot --filter-vcf $f_exac\n";
  print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot\n";
	print MAF "fi\n";
	print MAF "else\n";
    print MAF "     ".$run_script_path."remove_svtype.pl \${F_VCF_0} \${F_VCF_1}\n";
    print MAF "cd \${RUNDIR}\n";
    print MAF ". $script_dir/set_envvars\n";
    print MAF "     ".$run_script_path."vep_annotator_v1.1.pl ./vep.merged.input >&./vep.merged.log\n";
    print MAF "rm \${F_VCF_2}\n";
    print MAF "rm \${F_VEP_2}\n";
    print MAF "ln -s \${F_VCF_1} \${F_VCF_2}\n";
    print MAF "ln -s \${F_VEP_1} \${F_VEP_2}\n";
#    print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot --filter-vcf $f_exac\n"; 
	print MAF "     ".$run_script_path."vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $f_ref_annot\n";

	print MAF "fi\n";
	close MAF;
 	my $sh_file=$job_files_dir."/".$current_job_file;
    if($q_name eq "research-hpc")
    {
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>40000] rusage[mem=40000]\" -M 40000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";     }
    else {        $bsub_com = "bsub -q $q_name -n 1 -R \"select[mem>40000] rusage[mem=40000]\" -M 40000000 -w \"$hold_job_file\" -o $lsf_out -e $lsf_err bash $sh_file\n";   }
    print $bsub_com;
    system ($bsub_com);

    }
