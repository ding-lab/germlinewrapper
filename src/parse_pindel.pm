
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
