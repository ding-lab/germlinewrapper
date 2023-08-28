#!/usr/bin/env python

#title      :apply_AD_filter.py
#author      :Fernanda Martins Rodrigues (@fernanda); modified from Kuan Huang's script: https://github.com/ding-lab/PanCanAtlasGermline/blob/86b646f487f16c838a90b943764ec6b409dd47ed/analysis/process_files/germline/local/filter_VCF_AD.py
#date      :20191203

"""
This script takes as input a VCF file and filters variants according to a given allelic depth (AD) threshold.
    This script outputs a file with the suffix '.AD.<AD_threshold>.vcf'

    Usage:
            python apply_AD_filter.py [-h] -i <input VCF file> -AD <AD threshold> -O <output directory>

    Arguments:
            -i, --inputVCF:         input VCF file; gzip compressed
            -AD, --ADthreshold:     desired AD threshold; any variants with AD for the variant allele smaller than this threshold will be excluded
            -O, --outputDirectory:  directory to write output files to
            -h,  --help:           prints usage documentation

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
    parser.add_argument("-AD", "--ADthreshold", required=True, help="desired AD threshold; any variants with AD for the variant allele smaller than this threshold will be excluded")
    parser.add_argument("-o", "--outputfile", required=True, help="output file")

    args = vars(parser.parse_args())
    inputVCF = args["inputVCF"]
    ADthreshold = args["ADthreshold"]
    outputfile = args["outputfile"]

    return inputVCF, ADthreshold, outputfile


###############
## MAIN CODE ##
###############

def main():
    inputVCF, ADthreshold, outFile = argument_parser()
    AD_thres = float(ADthreshold)

    try:
        vcfF = gzip.open(inputVCF,"rt") 
    except IOError:
        print("VCF file does not exist!")

   # outFile_suffix = ".AD."+str(int(AD_thres))+".vcf"
   #inputFile_basename = inputVCF.split('/')[-1]
   #outFile = outputDirectory+inputFile_basename.replace(".vcf.gz",outFile_suffix)
    outF = open(outFile, "w")

    all_var = 0 # will count total number of variants in input file
    nonpass_AD_var = 0 # will count number of variants not passing the AD filter; i.e. variants with AD for the alternative allele smaller than the give threshold
    pass_var = 0 # will count number of variants passing AD filter

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
            info_var = str(var[7]).split(";") # get info field
            format_var = str(var[8]).split(":") # get format field
            geno_var = str(var[9]).split(":") # get genotype and AD fields
            AD_index = -1

            ### reference filter; if reference allele is "N", variant will be filtered out
            nonpass_ref = False
            if ref == "N":
                nonpass_ref = True

            ### AD filter (only keep variants with at least the AD threshold for the alternative allele)
            nonpass_AD = False
            for i in range(0,len(format_var)):
                if str(format_var[i]) == "AD":
                    AD_index = i

            AD = str(geno_var[AD_index])
            # for GATK and Pindel calls, we need to look at the second integer for the alternative allele
            if "," in AD:
                ref_alt_AD = AD.split(",")
                if int(ref_alt_AD[1]) < AD_thres:
                    nonpass_AD = True
            elif "/" in AD:
                nonpass_AD = True
            elif "." == AD:
                nonpass_AD = True
            # for VarScan calls
            else:
                if int(AD) < AD_thres:
                    nonpass_AD = True

            if not nonpass_ref and not nonpass_AD:
                outF.write(line + "\n")
                pass_var += 1
            else:
                nonpass_AD_var += 1

    # Filter summary:
    print("Number of total variants:", all_var)
    print("Number of variants failing AD filter of", AD_thres,":", nonpass_AD_var)
    print("Number of variants passing all thresholds:", pass_var)

    outF.close()


if __name__ == "__main__":
    main()

## END ##################################



