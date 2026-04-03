#!/usr/bin/env Rscript --vanilla

## title:            prep4manualreview_PECGS.R
## description:      Prepare CharGer outputs for Manual Review
## contact:          Fernanda Martins Rodrigues (fernanda@wustl.edu; mrodrigues.fernanda@gmail.com)
## last updated:     03/18/2026
## R version:        R version 4.5.3
## ggplot2:          ggplot2_4.0.2 
## dplyr:            dplyr_1.2.0
## tibble:           tibble_3.3.1

#### Import necessary libraries ====
## make sure you have these libraries installed before running (recommend using rocker/tidyverse:4 docker image)
library(optparse)
library(ggplot2)
library(tibble)
library(dplyr)


#### Setup user options ====
option_list = list(
  make_option(c("-w", "--work_dir"),
              type="character",
              default=NULL,
              help="working directory)",
              metavar="character"),
  make_option(c("-d", "--batch_date"),
              type="character",
              default=NULL,
              help="batch date in MMDDYYYY format",
              metavar="character"),
  make_option(c("-r", "--report_type"),
              type="character",
              default=NULL,
              help="N for normal only and TN for tumor-normal pairs",
              metavar="character"),
  make_option(c("-c", "--charger"),
              type="character",
              default=NULL,
              help="combined charger output - not filtered by frequency, with readcounts",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Check whether required arguments are provided
if (is.null(opt$work_dir)){
  print_help(opt_parser)
  stop("Path to working directory is required is required (--work_dir).n", call.=FALSE)
}

if (is.null(opt$batch_date)){
  print_help(opt_parser)
  stop("Batch date in MMDDYYYY format is required (--batch_date).n", call.=FALSE)
}

if (is.null(opt$report_type)){
  print_help(opt_parser)
  stop("Report type (N or TN) is required (--report_type).n", call.=FALSE)
}

if (is.null(opt$charger)){
  print_help(opt_parser)
  stop("Combined charger output - not filtered by frequency, with readcounts.n", call.=FALSE)
}


# set this to prevent memory issues
options(future.globals.maxSize= +Inf)

# set seed
set.seed(12345)

options(stringsAsFactors = F)

#### Main ====

## PE-CGS GERMLINE REVIEW

# Read in options

work_dir = opt$work_dir
batch_date = opt$batch_date
report_type = opt$report_type
charger = opt$charger

# Set working and outputs directory
dir.create(work_dir, recursive = T)
setwd(work_dir)

# Import CharGer results - not filtered by allele frequency
charged = read.delim(charger, sep = "\t", header = T)

# Import curated PECGS gene lists
genes = read.delim("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview/GeneLists/germline_lists_v11152023/PECGS_germline_master_gene_list_v3.tsv", header = T)

acmg_genes = genes[genes$ACMG == 1,]
pancan_genes = genes[genes$PanCancer == 1,]
mm_genes = genes[genes$MM == 1,]
chol_genes = genes[genes$CHOL == 1,]
crc_genes = genes[genes$CRC == 1,]

# Classify genes into pancancer, acmg, mm-specific, chol-specific, or crc-specific
charged$ACMG_SF_gene[ charged$HUGO_Symbol %in% acmg_genes$HGNC_Gene_Symbol ] = "yes"
charged$ACMG_SF_gene[ !charged$HUGO_Symbol %in% acmg_genes$HGNC_Gene_Symbol ] = "no"

charged$PanCancer_gene[ charged$HUGO_Symbol %in% pancan_genes$HGNC_Gene_Symbol ] = "yes"
charged$PanCancer_gene[ !charged$HUGO_Symbol %in% pancan_genes$HGNC_Gene_Symbol ] = "no"

charged$MM_specific_gene[ charged$HUGO_Symbol %in% mm_genes$HGNC_Gene_Symbol ] = "yes"
charged$MM_specific_gene[ !charged$HUGO_Symbol %in% mm_genes$HGNC_Gene_Symbol ] = "no"

charged$CHOL_specific_gene[ charged$HUGO_Symbol %in% chol_genes$HGNC_Gene_Symbol ] = "yes"
charged$CHOL_specific_gene[ !charged$HUGO_Symbol %in% chol_genes$HGNC_Gene_Symbol ] = "no"

charged$CRC_specific_gene[ charged$HUGO_Symbol %in% crc_genes$HGNC_Gene_Symbol ] = "yes"
charged$CRC_specific_gene[ !charged$HUGO_Symbol %in% crc_genes$HGNC_Gene_Symbol ] = "no"


# Filter out variants in genes not of interest
charged_cpgs = subset(charged, subset = (ACMG_SF_gene == "yes" |
                                           PanCancer_gene == "yes" |
                                           MM_specific_gene == "yes" |
                                           CHOL_specific_gene == "yes" |
                                           CRC_specific_gene == "yes"))

dim(charged_cpgs)
# [1] 42 43

# Mark rare and low frequency
charged_cpgs = add_column(charged_cpgs, .after = "Allele_Frequency", Allele_Frequency_class = NA)
charged_cpgs$Allele_Frequency_class[ charged_cpgs$Allele_Frequency <= 0.0005 ] = "Rare (MAF <= 0.05%)"
charged_cpgs$Allele_Frequency_class[ charged_cpgs$Allele_Frequency > 0.0005 & charged_cpgs$Allele_Frequency < 0.01 ] = "Low frequency (0.05% < MAF < 1%)"

# Filter out common variants
charged_cpgs = charged_cpgs[!is.na(charged_cpgs$Allele_Frequency_class),]
dim(charged_cpgs)
# [1] 39   44

## Add column to table, indicating which genes are oncogenes and with are tumor suppressors

tsgs.gsea = readLines("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview/GeneLists/tumor_suppressors_GSEA.txt")
oncogenes.gsea = readLines("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview/GeneLists/oncogenes_GSEA.txt")
volgestein = read.delim("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview/GeneLists/Volgestin2013Science_125genes_class.txt", header = TRUE)
cancer_gene_census = read.delim("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview/GeneLists/CancerGeneCensus_all_06022021.tsv", header=TRUE)

for (gene in charged_cpgs$HUGO_Symbol){
  if (gene %in% tsgs.gsea){
    charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Tumor suppressor gene"
  } else {
    if (gene %in% oncogenes.gsea){
      charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Oncogene"
    } else {
      if (gene %in% volgestein$Gene.Symbol){
        if (volgestein[volgestein$Gene.Symbol==gene,"Classification"]=="TSG"){
          charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Tumor suppressor gene"
        } else {
          if (volgestein[volgestein$Gene.Symbol==gene,"Classification"]=="Oncogene"){
            charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Oncogene"
          }
        }
      } else {
        if (gene %in% cancer_gene_census$Gene.Symbol){
          if (cancer_gene_census[cancer_gene_census$Gene.Symbol==gene, "Role.in.Cancer"] == "oncogene" | cancer_gene_census[cancer_gene_census$Gene.Symbol==gene, "Role.in.Cancer"] == "oncogene, fusion"){
            charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Oncogene"
          } else {
            if (cancer_gene_census[cancer_gene_census$Gene.Symbol==gene, "Role.in.Cancer"] == "TSG" | cancer_gene_census[cancer_gene_census$Gene.Symbol==gene, "Role.in.Cancer"] == "TSG, fusion"){
              charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Tumor suppressor gene"
            } else {
              if (cancer_gene_census[cancer_gene_census$Gene.Symbol==gene, "Role.in.Cancer"] == "oncogene, TSG" | cancer_gene_census[cancer_gene_census$Gene.Symbol==gene, "Role.in.Cancer"] == "oncogene, TSG, fusion"){
                charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Oncogene, Tumor suppressor gene"
              }
            }
          }
        } else {
          charged_cpgs[charged_cpgs$HUGO_Symbol==gene, "Gene_type"] = "Not classified"
        }
      }
    }
  }
}

## remove VCF_Details column to make file smaller:
charged_cpgs = subset(charged_cpgs, select = -c(VCF_Details) )

## PVS1 should only be applied to loss-of-function variants in genes where loss-of-function is a mechanism of disease (tumor suppressor genes in cancer, for example)
## therefore it should not be applied to oncogenes - this is a small bug in charger, I just fix it ad-hoc here ...
vars2fixPVS1 = unique(subset(charged_cpgs, subset=(grepl("PVS1", Positive_Evidence) == TRUE & grepl("Oncogene", Gene_type) == TRUE))$HGVSg)

# recalculate scores
for (var in vars2fixPVS1){
  row = which(charged_cpgs$HGVSg == var)
  charged_cpgs[row, "Positive_Evidence"] = gsub(",PVS1","",charged_cpgs[row, "Positive_Evidence"])
  charged_cpgs[row, "Positive_CharGer_Score"] = as.numeric(charged_cpgs[row, "Positive_CharGer_Score"])-8
  charged_cpgs[row, "CharGer_Score"] = as.numeric(charged_cpgs[row, "CharGer_Score"])-8
  score = charged_cpgs[row, "CharGer_Score"]
  if (score < 5){
    charged_cpgs[row, "CharGer_Classification"] = "Uncertain Signifcance"
  } else {
    if (score >= 5 & score < 9){
      charged_cpgs[row, "CharGer_Classification"] = "Likely Pathogenic"
    } else {
      if (score >= 9){
        charged_cpgs[row, "CharGer_Classification"] = "Pathogenic"
      }
    }
  }
  charged_cpgs[row, "CharGer_Summary"] = strsplit(charged_cpgs[row, "CharGer_Summary"]," -- ")[[1]][2]
}

# Reclassify variants - this follows Huang et al. 2018 methods
charged_cpgs$Overall_Classification = "Uncertain Significance"
charged_cpgs$Overall_Classification[charged_cpgs$CharGer_Score > 4] = "Prioritized VUS"
charged_cpgs$Overall_Classification[charged_cpgs$CharGer_Score > 8] = "Likely Pathogenic"
charged_cpgs$Overall_Classification[grep("PS1", charged_cpgs$Positive_Evidence)] = "Pathogenic"
charged_cpgs$Overall_Classification[charged_cpgs$ClinVar_Pathogenicity=="Pathogenic"] = "Pathogenic" 

# Add columns in which to manually input manual review decision and notes later
charged_cpgs$Manual_review = NA
charged_cpgs$Notes = NA
charged_cpgs$Cancer_related = NA

## Filter out known false positives - this is a list of variants that have been 
## proven to be false positives or benign based on current ClinVar release.
## We filter them out here so that our manual review workload is reduced. 

fp = read.delim("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview/4.ManualReview/false_positives_v03162026.tsv",
                header = T, sep = "\t")

charged_cpgs = subset(charged_cpgs, subset=(!HGVSg %in% fp$HGVSg))

# Manual review based on IGV images
if (report_type == "N"){
  write.table(charged_cpgs, paste0("all.samples.NormalOnly.", batch_date, ".charged2vcf.filtered.RCfiltered.4manualreview.tsv"), quote = F, row.names = F, sep = "\t")
} else if (report_type == "TN"){
  write.table(charged_cpgs, paste0("all.samples.TN.", batch_date, ".charged2vcf.filtered.RCfiltered.4manualreview.tsv"), quote = F, row.names = F, sep = "\t")
} else {
  print("Invalid report type. Will not add it to filename")
  write.table(charged_cpgs, paste0("all.samples.", batch_date, ".charged2vcf.filtered.RCfiltered.4manualreview.tsv"), quote = F, row.names = F, sep = "\t")
}


sessionInfo()

