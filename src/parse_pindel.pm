
sub parse_pindel {
    my $sample_name = shift;
    my $sample_full_path = shift;
    my $job_files_dir = shift;
    my $bsub = shift;
    my $REF = shift;
    my $perl = shift;
    my $gvip_dir = shift;
    my $pindel2vcf = shift;  #  $pindel_dir/pindel2vcf
    my $pindel_config = shift;

    $current_job_file = "j4_parse_pindel_g_".$sample_name.".sh";

    my $pindel_results = "$sample_full_path/pindel/pindel_out";
    my $filter_results = "$sample_full_path/pindel/filter_out";
    system("mkdir -p $filter_results");

    my $outlist="$filter_results/pindel.out.filelist";
    my $pin_var_file="$filter_results/pindel.out.raw";

## Pindel Filter - below is input into pindel_filter.v0.5
# lines below are added to data from $pindel_config
    die "$pindel_config does not exist\n" unless (-f $pindel_config);

    my $out = "$filter_results/pindel_filter.input";
    print("Copying $pindel_config to $out and appending\n");
    system("cp $pindel_config $out");

    open(OUT, ">>$out") or die $!;
    print OUT <<"EOF";
pindel.filter.pindel2vcf = $pindel2vcf
pindel.filter.variants_file = $pin_var_file
pindel.filter.REF = $REF
pindel.filter.date = 000000
EOF

# 1. Pull out all reads from pindel raw output with the label ChrID
#   http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
# 2. run pindel_filter.  Confirm what the output is - not documented.

    my $outfn = "$job_files_dir/$current_job_file";
    print("Writing to $outfn\n");
    open(OUT, ">$outfn") or die $!;
    print OUT <<"EOF";
#!/bin/bash

echo Collecting results in $pindel_results
find $pindel_results -name \'*_D\' -o -name \'*_SI\' -o -name \'*_INV\' -o -name \'*_TD\'  > $outlist
list=\$(xargs -a  $outlist)
cat \$list | grep ChrID > $pin_var_file

echo Running pindel_filter.v0.5.pl
$perl $gvip_dir/pindel_filter.v0.5.pl $filter_results/pindel_filter.input
EOF

    close OUT;
    my $bsub_com = "$bsub < $job_files_dir/$current_job_file\n";
    print("Executing:\n $bsub_com \n");

    system ( $bsub_com );

}

1;
