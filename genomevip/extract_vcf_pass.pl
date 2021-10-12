#!/usr/bin/env perl
# ---------------------
# @name GenomeVIP utility script to extract pass/failed VCF calls
# @author R. Jay Mashl <rmashl@wustl.edu>
# @version 0.1: original
#
# @syntax extract_vcf_pass.pl  input_vcf_file  output_pass_vcf_file  output_fail_vcf_file
# ---------------------
use warnings;
use strict;

my ($myorig, $mypass, $myfail) = @ARGV;

# read filter lines
open(IN, "< $myorig");
open(PASS, "> $mypass");
open(FAIL, "> $myfail");

while(<IN>) {
    if( /^#/ ) {
        print PASS $_;
        print FAIL $_;
    } else {
        my @a = split /\t/;
        if( $a[6] eq "PASS" ) {
            print PASS $_;
        } else {
            print FAIL $_;
        }
    }
}
close(FAIL);
close(PASS);
close(IN);
