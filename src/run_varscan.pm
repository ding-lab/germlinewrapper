# TODO: allow varscan, samtools to be defined in configuration file

sub run_varscan{
    my $IN_bam_N = shift;
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $varscan_config_snp = shift;
    my $varscan_config_indel = shift;

    my $varscan="/usr/local/VarScan.v2.3.8.jar";
    my $samtools="/usr/local/bin/samtools";


    die "File not found: $varscan_config_indel\n" if (! -e $varscan_config_indel);
    # ignore comments in varscan_config and convert newlines to spaces, so that all arguments are in one line
    my $varscan_snp_args=`grep -v "^#" $varscan_config_snp | tr '\n' ' '`;
    my $varscan_indel_args=`grep -v "^#" $varscan_config_indel | tr '\n' ' '`;

    my $workdir="$sample_full_path/varscan";
    system("mkdir -p $workdir");

    # Create a list of BAM files for varscan to use
    my $bam_list="$workdir/bamfilelist.inp";
    open(OUT, ">$bam_list") or die $!;
    print OUT "$IN_bam_N\n";
    close OUT;

    my $current_job_file = "j2_varscan_g_".$sample_name.".sh";
    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;

    my $outsnp="$workdir/$sample_name"."raw.snp.vcf";
    my $logsnp="$workdir/$sample_name"."raw.snp.log";
    my $outindel="$workdir/$sample_name"."raw.indel.vcf";
    my $logindel="$workdir/$sample_name"."raw.indel.log";

    print OUT <<"EOF";
#!/bin/bash
export JAVA_OPTS=\"-Xms256m -Xmx512m\"

$samtools mpileup -q 1 -Q 13 -B -f $REF -b $bam_list | java \${JAVA_OPTS} -jar $varscan mpileup2snp - $varscan_config_snp > $outsnp 2> $logsnp
$samtools mpileup -q 1 -Q 13 -B -f $REF -b $bam_list | java \${JAVA_OPTS} -jar $varscan mpileup2indel - $varscan_config_indel > $outindel 2> $logindel
EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");
    system ( $bsub_com );
}

1;
