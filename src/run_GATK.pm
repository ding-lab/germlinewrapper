
# Confirm input normal BAM file exists and is not empty
sub run_GATK {
    my $NBAM = shift;
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $gatk = shift;
    my $picard = shift;
    my $status_rg = shift;

    my $current_job_file = "j1_gatk_g_".$sample_name.".sh";
    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");


    my $workdir="$sample_full_path/gatk";
    system("mkdir -p $workdir");

    # Intermediate and output files
    # rg files are generated if status_rg = 0.  They are then removed
    # Their names are created by replacing trailing .bam of $NBAM with .rg.bam
    $NBAM_rg="${NBAM%.bam}.rg.bam";

    # rawvcf is output of Haplotype Caller step.
    $rawvcf="$workdir/$sample_name.raw.vcf";
    # rawvcf is labeled and copied to gvipvcf
    $gvipvcf="$workdir/$sample_name.gvip.vcf";

    # These two below are outputs
    $snvvcf="$workdir/$sample_name.snv.gvip.vcf";
    $indelvcf="$workdir/$sample_name.indel.gvip.vcf";

    # Step 1 is GATK HaplotypeCaller.  Optionally, call `picard AddOrReplaceReadGroups` based on value of status_rg
    # Input: $NBAM
    # Output: $rawvcf
    # if doing AddOrReplaceReadGroups,
    #   * Create NBAM_rg
    #   * Index NBAM_rg with samtools
    #   * run haplotype caller on NBAM_rg
    #   * Delete NBAM_rg and the index file NBAM_rg.bai

    # Step 2 labels the output from step 1 with GenomeVIP-specific annotation

    # Step 3 pulls out SNV and indel variants

    # Args to picard AddOrReplaceReadGroups.  This should be specified elsewhere
    my $picard_args = "RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20";

    # Args to HaplotypeCaller.  Thsi should be specified elsewhere
    my $haplotype_args = "-mbq 10 -rf DuplicateRead -rf UnmappedRead -stand_call_conf 10.0";

    my $step1;
    if ($status_rg) { # <srg> = bam having read group or not: 1, yes and 0, no (default 1)
        $step1 = <<'EOF'; 
            java  \${JAVA_OPTS} -jar $gatk -R $REF  -T HaplotypeCaller -I \${NBAM} -o \${rawvcf} $haplotype_args
EOF
    } else {
        $step1 = <<"EOF"; 
            java  \${JAVA_OPTS} -jar $picard AddOrReplaceReadGroups I=\${NBAM} O=\${NBAM_rg} $picard_args
            samtools index \${NBAM_rg}
            java  \${JAVA_OPTS} -jar $gatk -R $REF  -T HaplotypeCaller -I \${NBAM_rg} -o \${rawvcf} $haplotype_args
            rm \${NBAM_rg}
            rm \${NBAM_rg}.bai  # created by `samtools index`
EOF
    }

    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
        #!/bin/bash
        export JAVA_OPTS=\"-Xms256m -Xmx512m\"

        # Step 1
        $step1

        # Step 2
        $run_script_path/genomevip_label.pl GATK \${rawvcf} \${gvipvcf}

        # Step 3
        java \${JAVA_OPTS} -jar $gatk -R $REF -T SelectVariants  -V  $gvipvcf  -o $snvvcf  -selectType SNP -selectType MNP"."
        java \${JAVA_OPTS} -jar $gatk -R $REF -T SelectVariants  -V  $gvipvcf  -o $indelvcf  -selectType INDEL"."

        echo Written final result to $snvvcf and $indelvcf
EOF

    close OUT;

    $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");
    system ( $bsub_com );

}
1;
