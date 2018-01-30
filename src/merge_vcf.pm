
sub merge_vcf{
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $gatk = shift;

    $current_job_file = "j5_merge_vcf_g.".$sample_name.".sh";
    my $workdir = "$sample_full_path/merged/vcf";
    system("mkdir -p $workdir");

    my $gatk_snv_vcf="$sample_full_path/gatk/$sample_name.snv.gvip.filtered.vcf";
    my $gatk_indel_vcf="$sample_full_path/gatk/$sample_name.indel.gvip.filtered.vcf";
    my $varscan_snv_vcf="$sample_full_path/varscan/$sample_name.raw.snp.filtered.vcf";
    my $varscan_indel_vcf="$sample_full_path/varscan/$sample_name.raw.indel.filtered.vcf";
    my $pindel_vcf="$sample_full_path/pindel/filter_out/pindel.out.raw.CvgVafStrand_pass.Homopolymer_pass.vcf";
    my $merger_out="$workdir/merged.vcf";

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xmx10g\"
filter_gatk_varscan.pl \${RUNDIR} $sample_name
java \${JAVA_OPTS} -jar $gatk -R $REF -T CombineVariants -o \${merger_out} --variant:gsnp \${gatk_snv_vcf} --variant:gindel \${gatk_indel_vcf} --variant:vsnp \${varscan_snv_vcf} --variant:vindel \${varscan_indel_vcf} --variant:pindel \${pindel_vcf} -genotypeMergeOptions PRIORITIZE -priority gsnp,vsnp,gindel,vindel,pindel
EOF

    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}

1;
