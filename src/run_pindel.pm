
sub run_pindel{
    my $IN_bam_N = shift;
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $pindel = shift;  # path to pindel executable
    my $f_centromere = shift;

    $current_job_file = "j3_pindel_g_".$sample_name.".sh";

    my $workdir = "$sample_full_path/pindel/pindel_out";
    system("mkdir -p $workdir");

    my $config_fn = "$workdir/$sample_name.config";
    print("Writing to $config_fn\n");
    open(OUT, ">$config_fn") or die $!;
    print OUT <<"EOF";
$IN_bam_T\t500\t$sample_name.T  -- confirm this is deleted
$IN_bam_N\t500\t$sample_name.N
EOF

    my $pindel_args="-T 4 -m 6 -w 1"

    my $out = "$job_files_dir/$current_job_file";
    print("Writing to $out\n");
    open(OUT, ">$out") or die $!;
    print OUT <<"EOF";
#!/bin/bash

$pindel -f $REF -i $config_fn -o $workdir $pindel_args -J $f_centromere

EOF

    close OUT;

    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com );
}

1;
