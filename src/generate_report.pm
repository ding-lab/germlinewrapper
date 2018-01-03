

    print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
    $hold_job_file=$current_job_file;
    $current_job_file = "Run_report_gl_".$working_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    if(-e $lsf_out) 
	{
	`rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
	}
    open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
    print REPRUN "#!/bin/bash\n";
    print REPRUN "#BSUB -n 1\n";
    print REPRUN "#BSUB -R \"rusage[mem=40000]\"","\n";
    print REPRUN "#BSUB -M 40000000\n";
    #print REPRUN "#BSUB -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\'\n";
    #print STREKA "#BSUB -q long\n";
    #print REPRUN "#BSUB -q research-hpc\n";
    print REPRUN "#BSUB -q ding-lab\n";
    print REPRUN "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    print REPRUN "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    print REPRUN "#BSUB -J $current_job_file\n";
    print REPRUN "#BSUB -w \"$hold_job_file\"","\n";
    print REPRUN "      ".$run_script_path."generate_final_report.pl ".$run_dir."\n";
    close REPRUN;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

