#!/usr/bin/env python

# title   :get_site_lists.py
# author  :Fernanda Martins Rodrigues (@fernanda)
# date    :09-22-2019

"""
    This script takes a charger output file and outputs a file containing a list of regions that can be used as input for the tool bam-readcount (i.e. site lists file).

    This script will output a file with the suffix ".siteList".
    It expects charger output to have the first column as sample IDs, as it will output a siteList file for each sample in charger file.

    Usage:
            python get_site_lists.py [-h] -i <input VCF file> -O <output directory>

    Arguments:
            -i, --inputVCF:          input VCF file; gzip compressed
            -O, --outputDirectory:  directory to write output files to
            -h,  --help:              prints usage documentation

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
    parser.add_argument("-i", "--inputTSV", required=True, help=".tsv file output from charger; first column must be sample id (you might have to modify your charger output file prior to running this script")
    parser.add_argument("-o", "--out_dir", default=os.getcwd(), help="directory where the user would like the output file to be written to.")

    args = vars(parser.parse_args())
    input_tsv = args['inputTSV']
    out_dir = args['out_dir']

    if out_dir[-1] != '/':
        out_dir = out_dir + '/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    return input_tsv, out_dir

###############
## MAIN CODE ##
###############

def main():
    input_tsv, out_dir = argument_parser()

    # open charger file
    try:
        tsv = open(input_tsv,"r")
    except IOError:
        print("TSV file does not exist!")

    
    header=tsv.readline().strip().split('\t')
    if header[0] == "HUGO_Symbol":
        raise ValueError("Script expects first column of charger.tsv file to be SampleID. Your file's first column is HUGO_Symbol. Please edit your file and try again")
    else:
        lines=tsv.readlines()
        samples=[]

        for line in lines:
            line=line.strip().split('\t')
            samples.append(line[0])

        samples=set(samples)

        for sample in samples:
            outF=open(out_dir+sample+".siteList","w")
            for line in lines:
                line=line.strip().split('\t')
                if line[0] == sample:
                    outF.write(line[3]+'\t'+line[4]+'\t'+line[4]+'\n')
            outF.close()    


if __name__ == "__main__":
        main()

## END ##################################


