#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($dir, $prefix);

GetOptions(
    "d|dir=s"    => \$dir,
    "p|prefix=s" => \$prefix
) or die "Usage: perl add_rc.pl -d <dir> -p <prefix>\n";

die "Usage: perl add_rc.pl -d <dir> -p <prefix>\n" unless $dir && $prefix;

my $normal_file = "$dir/$prefix.N.rc.vaf";
my $tumor_file  = "$dir/$prefix.T.rc.vaf";

sub load_rc_vaf {
    my ($file) = @_;
    my %data;

    open(my $fh, "<", $file) or die "Cannot open $file: $!\n";
    while (<$fh>) {
        chomp;
        next if /^\s*$/;

        my @f = split /\t/;
        next unless @f >= 10;

        my ($sample, $chr, $start, $end, $ref, $alt, $total, $refc, $altc, $vaf) = @f[0..9];

        next if $sample eq "Sample";

        $chr =~ s/^chr//i;

        my $key = join("\t", $chr, $start, $ref, $alt);
        $data{$key} = [$refc, $altc, $vaf];
    }
    close $fh;

    return %data;
}

my %normal_rc = load_rc_vaf($normal_file);
my %tumor_rc  = load_rc_vaf($tumor_file);

my @inputs = (
    "$dir/$prefix.charged2vcf.filtered.tsv",
    "$dir/$prefix.charged2vcf.filtered.af0.0005.tsv"
);

foreach my $input_file (@inputs) {

    next unless -e $input_file;

    my $output_file = $input_file;
    $output_file =~ s/\.tsv$/.withrc.tsv/;

    my $filtered_output_file = $input_file;
    $filtered_output_file =~ s/\.tsv$/.withrc.filtered.tsv/;

    open(my $in,   "<", $input_file)           or die "Cannot open $input_file: $!\n";
    open(my $out,  ">", $output_file)          or die "Cannot write $output_file: $!\n";
    open(my $fout, ">", $filtered_output_file) or die "Cannot write $filtered_output_file: $!\n";

    my $line_num = 0;
    while (<$in>) {
        chomp;
        $line_num++;
        my @f = split /\t/;

        if ($line_num == 1) {
            my $header = $_ . "\tnref\tnalt\tnvaf\ttref\ttalt\ttvaf";
            print $out  $header, "\n";
            print $fout $header, "\n";
            next;
        }

        next unless @f >= 8;

        my ($chr, $pos, $ref, $alt) = ($f[3], $f[4], $f[6], $f[7]);
        $chr =~ s/^chr//i;

        my $key = join("\t", $chr, $pos, $ref, $alt);

        my ($nref, $nalt, $nvaf) = exists $normal_rc{$key} ? @{$normal_rc{$key}} : (".", ".", ".");
        my ($tref, $talt, $tvaf) = exists $tumor_rc{$key}  ? @{$tumor_rc{$key}}  : (".", ".", ".");

        my $out_line = join("\t", @f, $nref, $nalt, $nvaf, $tref, $talt, $tvaf);
        print $out $out_line, "\n";

        if ($nalt ne "." && $talt ne "." && $nvaf ne "." && $tvaf ne "." &&
            $nalt >= 5 && $talt >= 5 && $nvaf >= 0.2 && $tvaf >= 0.2) {
            print $fout $out_line, "\n";
        }
    }

    close $in;
    close $out;
    close $fout;

    print "Output written to: $output_file\n";
    print "Filtered output written to: $filtered_output_file\n";
}