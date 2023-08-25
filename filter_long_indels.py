#!/usr/bin/env python

#title   :filter_long_indels.py
#author  :Fernanda Martins Rodrigues (@fernanda);
#date    :20191908

"""
This script takes as input a VCF file and filters indels which ref and/or alt alleles are more than 100bp long.
        This step is good to remove misscalls from pindel.
        This script outputs a file with the suffix '.noLongIndels.vcf'

        Usage:
                python filter_long_indels.py [-h] -i <input VCF file> -O <output directory>

        Arguments:
                -i, --inputVCF:          input VCF file; gzip compressed
                -O, --outputDirectory:   directory to write output files to
                -h, --help:              prints usage documentation

"""


import sys
import argparse
import getopt
import gzip
import os

def argument_parser():
    # create parser
    parser = argparse.ArgumentParser(description=__doc__)
    # add arguments
    parser.add_argument("-i", "--inputVCF", required=True, help="input VCF file; gzip compressed")
    parser.add_argument("-O", "--outputDirectory", default=os.getcwd(), help="directory to write output files to")

    args = vars(parser.parse_args())
    inputVCF = args["inputVCF"]
    outputDirectory = args["outputDirectory"]

    if outputDirectory[-1] != '/':
        outputDirectory = outputDirectory + '/'

    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)

    return inputVCF, outputDirectory



###############
## MAIN CODE ##
###############


def main():
    vcfFile, outDir = argument_parser()

    try:
        vcfF = gzip.open(vcfFile,"rt")
    except IOError:
        print("VCF file does not exist!")

    outFile_suffix = ".noLongIndels.vcf"
    inputFile_basename = vcfFile.split('/')[-1]
    outFile = outDir+inputFile_basename.replace(".vcf.gz",outFile_suffix)
    outF = open(outFile, "w")

    all_var = 0 # will count total number of variants in input file
    nonpass_var = 0
    pass_var = 0

    # Start parsing input VCF file
    for line in vcfF:
        line=line.strip()
        # print the info lines to output file
        if line.startswith('#'):
            outF.write(line + "\n")
        else:
            var = line.split("\t")
            all_var += 1

            ref = str(var[3]) # get reference allele
            alt = str(var[4]) # get alternative allele

            ### if ref and/or alternative allele are bigger than 100 bp, variant will be filtered out.
            if len(ref) > 100:
                nonpass_var +=1
            elif len(alt) > 100:
                nonpass_var += 1
            else:
                pass_var += 1
                outF.write(line + "\n")

    # Filter summary:
    print("Number of total variants:", all_var)
    print("Number of variants failing indel filter:", nonpass_var)
    print("Number of variants passing indel filter:", pass_var)

    outF.close()


if __name__ == "__main__":
        main()

## END ##################################





