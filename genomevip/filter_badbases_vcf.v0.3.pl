#!/usr/bin/env perl
#----------------------------------
# @name Split out VCF records with IUPAC bases to separate file
# @author R. Jay Mashl <rmashl@wustl.edu>
# @version 0.3: do not include N afterall
# @version 0.2: allow base N for compatibility; append to badfile; allow optional explicit badfile name
# @version 0.1: original
#----------------------------------
use strict;
use warnings;

if( scalar @ARGV < 1 || scalar @ARGV > 2) {
    my $me = `basename $0`;
    chomp $me;
    print "\n     Syntax:  $me   <VcfFile>  [ <BadlinesOutfile> ]\n\n";
    exit 1;
}

my $vcf_file = $ARGV[0];
my $vcfstem   = `basename $vcf_file .vcf`;
chomp $vcfstem;

my $badfile   = "$vcfstem.iupac_ambiguous.vcf";
if( scalar @ARGV eq 2 ) {
    $badfile = $ARGV[1];
}


my @goodlines=();
my @badlines=();
open (VCF,   "<", $vcf_file  ) || die "Error: cannot open vcf file $vcf_file";
while( <VCF> ){
    my @a=split/\t/;
    if( /^#/ || ( $a[3] !~ tr/ACGT//c  &&  $a[4] !~ tr/ACGT,//c)) {
        push @goodlines, $_;
    } else {
        push @badlines, $_;
    }
}
close VCF;

if( scalar @badlines ) {
    print "Vcf file $vcf_file has records with ambiguous nucleotides. Putting these records into file $badfile\n";
    open( BAD, ">", $badfile ) || die "Error: cannot open/append file $badfile for output";
    foreach( @badlines ) {
        print BAD;
    }
    close BAD;
    print "Backing up $vcf_file to $vcf_file.old\n";
    system("cat $vcf_file > $vcf_file.old");
    open (CLEAN, ">", $vcf_file ) || die "Error: cannot open $vcf_file for rewriting";
    foreach( @goodlines ) {
        print CLEAN;
    }
    close CLEAN;
    print "Done.\n";
} else {
    print "Vcf file $vcf_file does not contain records with ambiguous nucleotides!\n";
}
