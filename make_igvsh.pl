#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $run_dir = "";

GetOptions(
    "rdir=s" => \$run_dir,
) or die "Usage: perl make_igvsh.pl --rdir <run_dir>\n";

die "Usage: perl make_igvsh.pl --rdir <run_dir>\n" unless $run_dir;

$run_dir =~ s/\/$//;

my $output_dir = $run_dir . ".manualreview";
my $input_tsv  = "$run_dir/all.samples.TN.htangl.charged2vcf.filtered.RCfiltered.4manualreview.tsv";

system("mkdir -p '$output_dir'") == 0 or die "Cannot create $output_dir\n";
system("mkdir -p '$output_dir/Snapshots_IGV'") == 0 or die "Cannot create $output_dir/Snapshots_IGV\n";

open(my $IN, "<", $input_tsv) or die "Cannot open $input_tsv: $!\n";

my $header = <$IN>;
die "Empty input file: $input_tsv\n" unless defined $header;
chomp $header;

my @cols = split /\t/, $header;
my %idx;
for my $i (0 .. $#cols) {
    $idx{$cols[$i]} = $i;
}

for my $req (qw/Sample Chromosome Start Stop/) {
    die "Missing required column: $req\n" unless exists $idx{$req};
}

my %sample_pos;

while (<$IN>) {
    chomp;
    next if /^\s*$/;

    my @f = split /\t/;

    my $sample = $f[$idx{'Sample'}];
    my $chr    = $f[$idx{'Chromosome'}];
    my $start  = $f[$idx{'Start'}];
    my $stop   = $f[$idx{'Stop'}];

    next unless defined $sample && defined $chr && defined $start && defined $stop;
    next if $sample eq "" || $chr eq "" || $start eq "" || $stop eq "";

    $chr =~ s/^chr//i;
    $chr = "chr" . $chr;

    $sample_pos{$sample}{$chr}{$start}{$stop} = 1;
}
close $IN;

sub chr_sort_key {
    my ($chr) = @_;
    $chr =~ s/^chr//i;
    return 23 if $chr eq 'X';
    return 24 if $chr eq 'Y';
    return 25 if $chr eq 'M' || $chr eq 'MT';
    return $chr if $chr =~ /^\d+$/;
    return 1000;
}

my $pad = 50;  # +/- 100 bp window

for my $sample (sort keys %sample_pos) {
    my $tumor_bam  = "$run_dir/$sample/$sample.T.bam";
    my $normal_bam = "$run_dir/$sample/$sample.N.bam";

    my $sh_file = "$output_dir/$sample.igv.sh";
    my $snapdir = "$output_dir/Snapshots_IGV/$sample";

    system("mkdir -p '$snapdir'") == 0 or die "Cannot create $snapdir\n";

    open(my $OUT, ">", $sh_file) or die "Cannot write $sh_file: $!\n";

    print $OUT "new\n";
    print $OUT "genome hg38\n";
    print $OUT "load $normal_bam\n\n";
    print $OUT "load $tumor_bam\n\n";

    for my $chr (
        sort {
            chr_sort_key($a) <=> chr_sort_key($b) || $a cmp $b
        } keys %{ $sample_pos{$sample} }
    ) {
        for my $start (sort { $a <=> $b } keys %{ $sample_pos{$sample}{$chr} }) {
            for my $stop (sort { $a <=> $b } keys %{ $sample_pos{$sample}{$chr}{$start} }) {

                # expand window
                my $left  = $start - $pad;
                my $right = $stop + $pad;
                $left = 1 if $left < 1;

                # custom snapshot name (no commas)
                my $png = "${chr}_${start}_${stop}.png";

                print $OUT "goto $chr $left $right\n";
                print $OUT "sort base\n";
                print $OUT "collapse\n";
                print $OUT "maxPanelHeight 2000\n";
                print $OUT "colorBy READ_STRAND\n";
                print $OUT "snapshotDirectory $snapdir/\n";
                print $OUT "snapshot $png\n";
            }
        }
    }

    close $OUT;
    print "Wrote $sh_file\n";
}