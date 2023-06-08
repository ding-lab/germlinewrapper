#!/usr/bin/env python

"""
This script takes as input an output vcf from germline wrapper and parses its contents in order to obtain one column per sample.
    The germline wrapper vcf output shows 3 columns per sample, one for each caller: gatk, pindel, varscan.
    This script parses input vcf according to the following rules:
    ### snv: union call from varscan and gatk.
    ### indel: callings from pindel or 2/3 callers; variants called by only gatk or only varscan will be filtered out.
    ### caller priority order: gsnp,vsnp,gindel,vindel,pindel

    This script outputs a file named merged.singleCaller.vcf

    Usage: python pickCaller.py -i <input vcf> -o <output file basename> -O <output directory>

"""
import sys
import argparse
import os

### Define arguments ###

def argument_parser():
    # create parser
    parser = argparse.ArgumentParser(description=__doc__)
    # add arguments
    parser.add_argument("-i", "--input_vcf", required=True, help="input vcf file; output from germline wrapper")
    parser.add_argument("-o", "--out_basename", required=True, help="basename for output vcf file")
    parser.add_argument("-O", "--out_dir", default=os.getcwd(), help="output directory; optional")

    args = vars(parser.parse_args())

    input_vcf = args['input_vcf']
    out_basename = args['out_basename']
    out_dir = args['out_dir']

    if out_dir[-1] != '/':
        out_dir = out_dir + '/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    return input_vcf, out_basename, out_dir



### MAIN CODE ###

def main():
    input_vcf, out_basename, out_dir = argument_parser()

    vcf = open(input_vcf)

    # open new vcf file to write output to:
    new_vcf=open(out_dir+out_basename+".singleCaller.vcf","w")

    # parse old vcf file
    for line in vcf:
        if line.startswith("##"):
            new_vcf.write(line)
        elif line.startswith("#CHROM"):
            line=line.strip().split("\t")
            new_vcf.write("\t".join(line[0:9])+"\t"+out_basename+"\n")
        else:
            line=line.strip().split("\t")
            genotypes=line[9:] # all three callers
            # remove empty genotypes
            genotypes=[g for g in genotypes if not "./." in g]

            # first parse snps (call from either gatk or varscan suffices; if both, take gatk over varscan)
            if "snp" in line[7]:
                info=line[7].strip().split(";")
                if info[-1] == "set=gsnp-vsnp":
                    info[-1]="set=snp"
                    info=";".join(info)
                    gLength1=len(genotypes[0])
                    gLength2=len(genotypes[1])
                    if gLength1 < gLength2: # varscan has more format fields than gatk, so string is longer
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[0]+"\n")
                    else:
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[1]+"\n")
                else:
                    info[-1]="set=snp"
                    info=";".join(info)
                    new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[0]+"\n")
            # now, parse indels (call from pindel or both gatk and varscan; if all three, take gatk over the others)
            elif "indel" in line[7]:
                info=line[7].strip().split(";")
                if info[-1] == "set=gindel-vindel-pindel":
                    info[-1]="set=indel"
                    info=";".join(info)
                    gLength1=len(genotypes[0])
                    gLength2=len(genotypes[1])
                    gLength3=len(genotypes[2])
                    if (gLength1>gLength2>gLength3) == True or (gLength1<gLength2<gLength3):
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[1]+"\n")
                    elif (gLength2>gLength1>gLength3) == True or (gLength2<gLength1<gLength3):
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[0]+"\n")
                    elif (gLength1>gLength3>gLength2) == True or (gLength1<gLength3<gLength2):
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[2]+"\n")
                elif info[-1] == "set=gindel-vindel":
                    info[-1]="set=indel"
                    info=";".join(info)
                    gLength1=len(genotypes[0])
                    gLength2=len(genotypes[1])
                    if gLength1 < gLength2: # varscan has more format fields than gatk, so string is longer
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[0]+"\n")
                    else:
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[1]+"\n")
                elif info[-1] == "set=vindel-pindel" or info[-1] == "set=gindel-pindel":
                    info[-1]="set=indel"
                    info=";".join(info)
                    gLength1=len(genotypes[0])
                    gLength2=len(genotypes[1])
                    if gLength1 < gLength2: # varscan has more format fields than gatk, so string is longer
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[1]+"\n")
                    else:
                        new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[0]+"\n")
                elif info[-1] == "set=pindel":
                    info[-1]="set=indel"
                    info=";".join(info)
                    new_vcf.write("\t".join(line[0:7])+"\t"+info+"\t"+line[8]+"\t"+genotypes[0]+"\n")
    
    vcf.close()
    new_vcf.close()

    return 0

if __name__=='__main__':
    main()

### END ###
