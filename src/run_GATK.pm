
sub bsub_gatk{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_gatk_g_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    if(-e $lsf_out) {
        `rm $lsf_out`;
        `rm $lsf_err`;
        `rm $current_job_file`;
    }
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
