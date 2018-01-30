
# From somaticwrapper/src/run_vep.pl
sub write_vep_input {
    my $config_fn = shift;
    my $module = shift;  # e.g. varscan.vep or strelka.vep
    my $vcf = shift;
    my $output = shift;
    my $vep_cmd = shift;
    my $cache_dir = shift;
    my $REF = shift;
    my $assembly = shift;
    my $use_vep_db = shift;  # 1 for testing/demo, 0 for production
    my $output_vep = shift;  # output annotated vep rather than vcf format after merge step.  add suffix 'vep' to output

    if ($output_vep) {
        $output = "$output.vep";
    }

    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$module.vcf = $vcf
$module.output = $output
$module.vep_cmd = $vep_cmd
$module.cachedir = $cache_dir
$module.reffasta = $REF
$module.assembly = $assembly
$module.usedb = $use_vep_db
$module.output_vep = $output_vep
EOF
}

sub bsub_vcf_2_maf{
    my $IN_bam_N = $sample_full_path."/".$sample_name.".N.bam";
    $vep
    $vep_cache
$REF
$assembly
$f_exac
  
    $current_job_file = "j6_vcf_2_maf.".$sample_name.".sh";

    my $workdir = "$sample_full_path/merged/maf";
    system("mkdir -p $workdir");


    my $F_VCF_0=$workdir."/merged.vcf";
    my $F_VCF_1=$workdir."/merged.1.vcf";
    my $F_VCF_2=$workdir."/".$sample_name.".vcf";
    my $F_VEP_1=$workdir."/merged.VEP.vcf";
    my $F_VEP_2=$workdir."/".$sample_name.".vep.vcf";
    my $F_maf=$workdir."/".$sample_name.".maf";

    my $vep_config = $workdir."/vep.merged.input";
    write_vep_input(
        "$workdir/vep.merged.input",
        "merged.vep",
        $F_VEP_0,
        $F_VEP_1,
        "$filter_results/merged.VEP.vcf",
        $vep_cmd, $cache_dir, $REF, $assembly, $use_vep_db, $output_vep);

#cat > $vep_config <<EOF
#merged.vep.vcf = $F_VCF_0
#merged.vep.output = $F_VEP_1
#merged.vep.vep_cmd = $vep
#merged.vep.cachedir = $vep_cache
#merged.vep.reffasta = $REF
#merged.vep.assembly = $assembly
#EOF


    open(MAF, ">$job_files_dir/$current_job_file") or die $!;
#!/bin/bash

# this appears to simply remove the string "SVTYPE=" from 0 and write to 1
remove_svtype.pl \${F_VCF_0} \${F_VCF_1}

#. /gscmnt/gc2525/dinglab/rmashl/Software/perl/set_envvars
    # Contents of above
    # export PERL_PATH=/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0
    # export PATH=$PERL_PATH/bin:$PATH
    # export PERL_BIN=$PERL_PATH/bin/perl
    # export PERL5LIB=$PERL_PATH/lib/perl5:$PERL_PATH/lib/perl5/x86_64-linux:$PERL5LIB

vep_annotator.pl ./vep.merged.input >& ./vep.merged.log
rm \${F_VCF_2}
rm \${F_VEP_2}
ln -s \${F_VCF_1} \${F_VCF_2}
ln -s \${F_VEP_1} \${F_VEP_2}
vcf2maf.pl --input-vcf \${F_VCF_2} --output-maf \${F_maf} --tumor-id $sample_name\_T --normal-id $sample_name\_N --ref-fasta $REF --filter-vcf $f_exac

# splice_site_check.pl \${F_maf}  # According to Song 1/11/18 no longer needed


    close MAF;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);

    }
