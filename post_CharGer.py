#!/usr/bin/env python

#title   :post_CharGer.py
#author	 :Fernanda Martins Rodrigues (fernanda@wustl.edu)
#date	 :12282023

"""
	This script fixes the HGVSc and HGVSp annotations output from CharGer, which can be problematic, and pulls allele frequency information from
	1000 Genomes, ExAC, and gnomAD.
	It pulls these annotations from the VEP annotation, present in the VCF_Details colum of the CharGer output.
	It also parses the CharGer output file and its corresponding VCF file, combining genotype information from VCF to CharGer output.

	Finally, this script will filter CharGer output as follows:
		- Filter out any variants classified by CharGer as Benign or Likely Benign;
		- Filter out synonymous variants;
		- Filter out variants with CharGer score < 4.

	Usage:
		python post_CharGer.py [-h] -i <input VCF file> -c <input CharGer tsv file> -d <disease-optional> -s <sample ID> -a <MAF threshold> -o <basename for output file> -O <output directory>

	Arguments:
		-i, --inputVCF:			input VCF file
		-c, --charger:			charger output tsv file corresponding to input VCF file
		-d, --disease:			disease; optional
		-s, --sampleID:	  		sample ID to add to output file and output file name
		-a, --maf:				MAF threshold
		-o, --outBasename:	  	output file basename
		-O, --outputDirectory:  directory to write output files to

"""

import sys
import argparse
import getopt
import os
import re
import gzip

def argument_parser():
	# create parser
	parser = argparse.ArgumentParser(description=__doc__)
	# add arguments
	parser.add_argument("-i", "--inputVCF", required=True, help="input VCF file")
	parser.add_argument("-c", "--charger", required=True, help="charger output tsv file corresponding to input VCF file")
	parser.add_argument("-d", "--disease", required=False, help="optional disease parameter; this will add a column 'Disease' to the output file")
	parser.add_argument("-s", "--sampleID", required=True, help="sample ID to add to output file and output file name")
	parser.add_argument("-a", "--maf", required=True, help="minor allele frequency (MAF) threshold to filter charger output")
	parser.add_argument("-o", "--outBasename", required=True, help="output file basename")
	parser.add_argument("-O", "--outputDirectory", required=True, help="directory to write output files to")

	args = vars(parser.parse_args())
	inputVCF = args["inputVCF"]
	charger = args["charger"]
	disease = args["disease"]
	sampleID = args["sampleID"]
	maf = args["maf"]
	outBasename = args["outBasename"]
	outputDirectory = args["outputDirectory"]

	if outputDirectory[-1] != '/':
		outputDirectory = outputDirectory + '/'

	if not os.path.exists(outputDirectory):
		os.makedirs(outputDirectory)

	return inputVCF, charger, disease, sampleID, maf, outBasename, outputDirectory


###############
## FUNCTIONS ##
###############

def get_hgvsc(charger_header, charger_line_list):
	'''
	Gathers correct HGVSc annotation from VCF_Details column in CharGer output.
	'''
	# get VCF info in charger output
	vcf_details = charger_line_list[charger_header.index("VCF_Details")].split("|")

	# get correct ENST (transcript ID) to search for correct HGVSc annotation in VCF info:
	hgvsc = charger_line_list[charger_header.index("HGVSc")]
	enst = hgvsc.split(":")[0]

	# find correct HGVSc:
	enst_r = re.compile(enst+".*:c.")
	try:
		hgvsc_vep = list(filter(enst_r.match, vcf_details))[0] # HGVSc annotation with ENST (transcript) ID
		hgvsc_short = hgvsc_vep.split(":")[1] # HGVSp annotation without transcript ID
	except:
		hgvsc_vep = enst+':c.'
		hgvsc_short = ':c.'

	return hgvsc_vep, hgvsc_short


def get_hgvsp(charger_header, charger_line_list, hgvsc_vep):
	'''
	Gathers correct HGVSp annotation from VCF_Details column in CharGer output.
	hgvsc_vep is output of get_hgvsc() function
	'''

	# define amino-acid dictionary
	aa_code = {'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe':
	'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L', 'Met':
	'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S', 'Thr':
	'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y', 'Ter': '*'}

	# get VCF info from charger output
	vcf_details = charger_line_list[charger_header.index("VCF_Details")].split("|")
	var_type = charger_line_list[charger_header.index("Variant_Classification")]

	# get ENSP (protein ID) to search for correct HGVSp in VCF information
	hgvsp = charger_line_list[charger_header.index("HGVSp")]
	ensp = hgvsp.split(":")[0]
	ensp_r = re.compile(ensp+".*:p.")

	try:
		hgvsp_vep = list(filter(ensp_r.match, vcf_details))[0]
	except:
		hgvsp_vep = ensp+":p."
		hgvsp_short = "p."
		return hgvsp_vep, hgvsp_short

	# convert 3 letter aminoacid to 1 letter in hgvsp annotation
	p_annot = hgvsp_vep.split(":")[1]
	p_annot_list = list(re.findall(r"(.*:*p\.)(\D+)(\d+)(\D+)",hgvsp_vep,re.IGNORECASE)[0])

	try:
		p_annot_list[1] = aa_code[ p_annot_list[1] ]
	except:
		hgvsp_vep = ensp+":p."
		hgvsp_short = "p."
		return hgvsp_vep, hgvsp_short

	if p_annot_list[3] in aa_code:
		p_annot_list[3] = aa_code[ p_annot_list[3] ]
	elif p_annot_list[3] == "?":
		hgvsp_pos = vcf_details[vcf_details.index(hgvsc_vep)+4]
		hgvsp_aa = vcf_details[vcf_details.index(hgvsc_vep)+5].split("/")
		if len(hgvsp_aa) == 1 and hgvsp_aa[0] == '':
			print("AA FIELD IS EMPTY FOR ", p_annot_list)
		else:
			p_annot_list[1] = hgvsp_aa[0]
			p_annot_list[3] = hgvsp_aa[1]
	elif var_type == "inframe_deletion":
		p_annot_list[3] = "del"
	elif var_type == "inframe_insertion":
		p_annot_list[3] = "ins"
	elif var_type == "synonymous_variant":
		p_annot_list[3] = p_annot_list[1]
	else:
		p_annot_list[3] = "fs"
	hgvsp_vep = "".join(p_annot_list)
	hgvsp_short = hgvsp_vep.split(":")[1]

	return hgvsp_vep, hgvsp_short


def get_MAFs(charger_header, vep_header, charger_line_list):
	'''
	Get MAF from 100 Genomes, ExAC, and gnomAD from VEP annotation in CharGer file (VCF_Details column)

	'''
	# get VCF info from charger output
	vcf_details = charger_line_list[charger_header.index("VCF_Details")].split("|")

	# define VEP info field information: (this is according to vep 100 using the --everything parameter)
	vep=vep_header.split('|')

	# parse VEP
	if len(vcf_details[vep.index("AF")]) == 0 or vcf_details[vep.index("AF")] == '':
		kgMAF = 'NA'
	else:
		kgMAF = vcf_details[vep.index("AF")]

	if len(vcf_details[vep.index("ExAC_Adj_AF")]) == 0 or vcf_details[vep.index("ExAC_Adj_AF")] == '':
		exacMAF = 'NA'
	else:
		exacMAF = vcf_details[vep.index("ExAC_Adj_AF")]

	if len(vcf_details[vep.index("gnomAD_AF")]) == 0 or vcf_details[vep.index("gnomAD_AF")] == '':
		gnomadMAF = 'NA'
	else:
		gnomadMAF = vcf_details[vep.index("gnomAD_AF")]


	return kgMAF, exacMAF, gnomadMAF


def get_charger_varcoord(charger_header, charger_line_list):
	'''
	Get variant coordinates from CharGer.
	This will be used to match the variant to corresponding VCF
	'''

	# get VCF info from charger output
	vcf_details = charger_line_list[charger_header.index("VCF_Details")]

	# get variant coordinates to match VCF
	if "," in vcf_details.split(':')[8][2:-2]:
		alt1=vcf_details.split(':')[8].split(",")[0][2:-1]
		alt2=vcf_details.split(':')[8].split(",")[1][1:-2]
		var_coord=[vcf_details.split(':')[0], vcf_details.split(':')[2], vcf_details.split(':')[6], alt1+','+alt2]
	else:
		var_coord=[vcf_details.split(':')[0], vcf_details.split(':')[2], vcf_details.split(':')[6], vcf_details.split(':')[8][2:-2]]

	return var_coord



def get_genotype(vcf_lines):
	'''
	Functions loops through lines in a given VCF and parses out genotype information
	'''
	vcf_vars=[]
	vcf_genotypes=[]

	for line in vcf_lines:
		if line.startswith( '#' ):
			continue
		else:
			line=line.strip().split('\t')
			vcf_vars.append([line[0], line[1], line[3], line[4]])
			vcf_genotypes.append(line[9])

	return vcf_vars, vcf_genotypes


def get_vep_header(vcf_lines):
	'''
	Get VEP header from VCF
	'''
	for line in vcf_lines:
		if line.startswith('#') and "from Ensembl VEP. Format:" in line:
			line = line.strip().split("from Ensembl VEP. Format: ")
			vep=line[1][0:-2]
			print(vep)

	return vep

###############
## MAIN CODE ##
###############

def main():
	inputVCF, charger, disease, sampleID, maf, outBasename, outputDirectory = argument_parser()

	# Open input files
	if ".gz" in inputVCF:
		try:
			vcf=gzip.open(inputVCF, 'rt')
		except IOError:
			print("Input vcf.gz file does not exist!")
			sys.exit(1)
	else:
		try:
			vcf=open(inputVCF, 'r')
		except IOError:
			print("Input vcf file does not exist!")
			sys.exit(1)
	try:
		tsv=open(charger)
	except IOError:
		print("Input CharGer file does not exist!")
		sys.exit(1)

	### 1. Map charger to vcf ###
	# Open output file
	outFilename1=outputDirectory+sampleID+".charged2vcf.tsv"
	charger2vcf=open(outFilename1, "w")

	# Header for output file
	charger_header=tsv.readline().strip().split("\t")
	if disease:
		out_header = "Sample\tDisease\tHUGO_Symbol\tChromosome\tStart\tStop\tReference\tAlternate\tVariant_Classification\tHGVSg\tHGVSc\tHGVSc_short\tHGVSp\tHGVSp_short\tAllele_Frequency\tVEP_Most_Severe_Consequence\tPositive_Evidence\tNegative_Evidence\tPositive_CharGer_Score\tNegative_CharGer_Score\tCharGer_Score\tClinVar_Pathogenicity\tACMG_Classification\tCharGer_Classification\tPubMed_Link\tClinVar_Traits\tCharGer_Summary\tVCF_Details\tGenotype\t1KGenomes_AF\tExAC_AF\tgnomAD_AF\n"
	else:
		out_header = "Sample\tHUGO_Symbol\tChromosome\tStart\tStop\tReference\tAlternate\tVariant_Classification\tHGVSg\tHGVSc\tHGVSc_short\tHGVSp\tHGVSp_short\tAllele_Frequency\tVEP_Most_Severe_Consequence\tPositive_Evidence\tNegative_Evidence\tPositive_CharGer_Score\tNegative_CharGer_Score\tCharGer_Score\tClinVar_Pathogenicity\tACMG_Classification\tCharGer_Classification\tPubMed_Link\tClinVar_Traits\tCharGer_Summary\tVCF_Details\tGenotype\t1KGenomes_AF\tExAC_AF\tgnomAD_AF\n"

	charger2vcf.write(out_header)


	# Parse VCF file and store genotype information for each variant
	vcf_lines = vcf.readlines()
	vcf_vars, vcf_genotypes = get_genotype(vcf_lines)
	vep = get_vep_header(vcf_lines) # also get vep for use later

	# Parse CharGer file and fix HVSc and HGVSp annotation

	charger_lines=tsv.readlines()

	for line in charger_lines:
		line=line.strip().split('\t')

		# find correct HGVSc annotation
		hgvsc_vep, hgvsc_short = get_hgvsc(charger_header, line)

		# find correct HGVSp annotation
		hgvsp_vep, hgvsp_short = get_hgvsp(charger_header, line, hgvsc_vep)

		# store variant coordinate information for matching to VCF
		var_charger = get_charger_varcoord(charger_header, line)

		# find genotype for variant already parsed from VCF
		if var_charger not in vcf_vars:
			print("Variant doesn't match VCF:", var_charger)
		for i in vcf_vars:
			if i == var_charger:
				genotype = vcf_genotypes[vcf_vars.index(i)]

		# get 1000 genomes, exac, and gnomad MAFs
		kgMAF, exacMAF, gnomadMAF = get_MAFs(charger_header, vep, line)

		# write out to output file
		if disease:
			out_line = sampleID + '\t' + disease + '\t' + '\t'.join(line[0:8])+'\t'+hgvsc_vep+'\t'+hgvsc_short+'\t'+hgvsp_vep+'\t'+hgvsp_short+'\t'+'\t'.join(line[10:])+'\t'+genotype+'\t'+kgMAF+'\t'+exacMAF+'\t'+gnomadMAF+'\n'
		else:
			out_line = sampleID + '\t' + '\t'.join(line[0:8])+'\t'+hgvsc_vep+'\t'+hgvsc_short+'\t'+hgvsp_vep+'\t'+hgvsp_short+'\t'+'\t'.join(line[10:])+'\t'+genotype+'\t'+kgMAF+'\t'+exacMAF+'\t'+gnomadMAF+'\n'

		charger2vcf.write(out_line)

	charger2vcf.close()

	### 2. Filter CharGer output ###
	# Open mapped file
	charger2vcf = outFilename1
	try:
		tsv=open(charger2vcf)
	except IOError:
		print("Input CharGer file does not exist!")

	# Open output file
	outFilename2 = outputDirectory+outBasename+".charged2vcf.filtered.tsv"
	outFilename2_AF = outputDirectory+outBasename+".charged2vcf.filtered.af"+str(maf)+".tsv"

	outFile=open(outFilename2, "w")
	outFileAF=open(outFilename2_AF, "w")

	# Parse and filter CharGer file
	charger_header = tsv.readline()
	outFile.write(charger_header)
	outFileAF.write(charger_header)

	header = charger_header.strip().split("\t")

	for line in tsv:
		info = line.strip().split("\t")
		# filter out variants classified by CharGer as Benign or Likely Benign:
		if "Benign" not in info[header.index("CharGer_Classification")]:
			# filter out synonymous variants and variants with CharGer score < 4
			if info[header.index("Variant_Classification")] != "synonymous_variant" and int(info[header.index("CharGer_Score")]) > 4:
				outFile.write(line)
				# filter out common variants
				if float(info[header.index("Allele_Frequency")]) <= float(maf):
					outFileAF.write(line)

	outFile.close()
	outFileAF.close()

if __name__ == "__main__":
	main()


### END ###



