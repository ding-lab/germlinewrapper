
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
    if(-e $lsf_out) {
        `rm $lsf_out`;
        `rm $lsf_err`;
        `rm $current_job_file`;
    }

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
