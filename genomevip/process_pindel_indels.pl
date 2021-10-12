#!/usr/bin/env perl
#----------------------------------
# @name Split out indels from pindel vcf output
# @author R. Jay Mashl
# @version 0.1
#----------------------------------
use strict;
use warnings;

use Getopt::Long;

my $vcf_file;
GetOptions (
    'variants=s'  =>  \$vcf_file
    ) || die "Error in command line arguments\n";

my $vcfstem   = `basename $vcf_file .vcf`;
chomp $vcfstem;
my $indelfile = "$vcfstem.indel.vcf";
my $otherfile = "$vcfstem.nonindel.vcf";

open (VCF,   "<", $vcf_file  ) || die "Error: cannot open vcf file $vcf_file";
open (INDEL, ">", $indelfile ) || die "Error: cannot open indel file $indelfile for output";
open (OTHER, ">", $otherfile ) || die "Error: cannot open nonindel file $otherfile for output";

while( <VCF> ){
    if( /^#/ ){
        print INDEL;
        print OTHER;
    } else {
	if( /SVTYPE=(INS|DEL)/ ){
	    print INDEL;
	} else {
	    print OTHER;
	}
    }
}
close OTHER;
close INDEL;
close VCF;
