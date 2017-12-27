
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
