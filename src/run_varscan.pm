
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
    if(-e $lsf_out) {
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
	print VARSCAN "ncols=\$(echo \"3*( \$(wc -l < \$BAMLIST) +1)\"|bc)\n";
    print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar mpileup2snp  -  --p-value  0.10   --min-coverage  3   --min-var-freq  0.08   --min-reads2  2   --min-avg-qual  15   --min-freq-for-hom  0.75   --strand-filter  1   --output-vcf  1   > \${outsnp}  2> \${logsnp}\n";   
 	print VARSCAN "\${SAMTOOLS_DIR}/samtools mpileup -q 1 -Q 13 -B -f $h37_REF -b \${BAMLIST} | awk -v ncols=\$ncols \'NF==ncols\' | java \${JAVA_OPTS} -jar \${VARSCAN_DIR}/VarScan.jar mpileup2indel  -  --p-value  0.10   --min-coverage  3   --min-var-freq  0.20   --min-reads2  2   --min-avg-qual  15   --min-freq-for-hom  0.75   --strand-filter  1   --output-vcf  1   > \${outindel}  2> \${logindel}\n";
    close VARSCAN;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

	}
