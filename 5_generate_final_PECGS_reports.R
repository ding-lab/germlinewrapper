#!/usr/bin/env Rscript --vanilla

## title:            5_generate_final_PECGS_reports.R
## description:      Prepare final PE-CGS germline reports
## contact:          Fernanda Martins Rodrigues (fernanda@wustl.edu; mrodrigues.fernanda@gmail.com)
## last updated:     03/19/2026
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
library(stringi)
library(tidyverse)
library(openxlsx)


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
  make_option(c("-s", "--samples2cancer"),
              type="character",
              default="all_processed_sample_ids_2cancer.tsv",
              help=".tsv file containing case_id\tsample_id\tdisease_type of all samples processed up to the current batch, including the current batch.",
              metavar="character"),
  make_option(c("-N", "--current_normalonly_cases"),
              type="character",
              default=NULL,
              help=".tsv file containing case ids for normal only samples for the CURRENT BATCH ONLY.",
              metavar="character"),
  make_option(c("-T", "--current_tumornormal_cases"),
              type="character",
              default=NULL,
              help=".tsv file containing case ids for tumor-normal samples for the CURRENT BATCH ONLY.",
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

# set this to prevent memory issues
options(future.globals.maxSize= +Inf)

# set seed
set.seed(12345)

options(stringsAsFactors = F)

#### Main ====

# Read in options

work_dir = opt$work_dir
batch_date = opt$batch_date
samples2cancer = opt$samples2cancer

# Set working and outputs directory
dir.create(work_dir, recursive = T)
dir.create(paste0(work_dir,"/5.CohortFrequencies/", batch_date), recursive = T)
#setwd(work_dir)

# Set case lists
if (is.null(opt$current_normalonly_cases) || opt$current_normalonly_cases == "") {
  print("NO NORMAL ONLY CASE LIST GIVEN, WILL LOOK FOR CURRENT BATCH FILE...")
  file_path = paste0("case_ids_", batch_date, "_GWprocessing_NormalOnly.tsv")
} else {
  file_path = opt$current_normalonly_cases
}

if (file.exists(file_path)) {
  current_normalonly_cases <- readLines(file_path)
} else {
  print("File not found. Proceeding with empty case list.")
  current_normalonly_cases <- c()
}
  

if (is.null(opt$current_tumornormal_cases) || opt$current_tumornormal_cases == "") {
  print("NO TN CASE LIST GIVEN, WILL LOOK FOR CURRENT BATCH FILE...")
  file_path = paste0("case_ids_", batch_date, "_GWprocessing_TN.tsv")
} else {
  file_path = opt$current_tumornormal_cases
}

if (file.exists(file_path)) {
  current_tumornormal_cases <- readLines(file_path)
} else {
  print("File not found. Proceeding with empty case list.")
  current_tumornormal_cases <- c()
}

current_cases = unique(c(current_normalonly_cases, current_tumornormal_cases))

#### Calculate cohort frequencies ====

# Import all ids processed up to the current batch
samples = read.delim(samples2cancer, header = T, sep = "\t")

# Generate master table of all ACMG, low freq, and rare vars
file_list = list.files(paste0(work_dir,"/4.ManualReview"), recursive = T, pattern = "manuallyReviewed.tsv", full.names = T)

all_vars = data.frame()

for (file in file_list){
  df = read.delim(file, header = T, sep = "\t")
  if ("Cancer_related" %in% colnames(df)){
    df = subset(df, select = -c(Cancer_related) )
  }
  all_vars = rbind(all_vars, df)
}

# Filter out variants failing manual review
all_vars = all_vars[all_vars$Manual_review == "PASS",]

# Add case ID
cases = unlist( lapply( strsplit(all_vars$Sample, "-"), "[", 1))
all_vars = add_column(all_vars, Case_ID = cases, .after = "Sample")

# Summarize results at the case level
all_vars = all_vars[,2:48]
all_vars = distinct(all_vars)
dim(all_vars)

# Take highest VAF
all_vars$mut_id = paste0(all_vars$Case_ID,"_",all_vars$HGVSg)
all_vars = all_vars %>% group_by(mut_id) %>% slice_max(tvaf, n = 1)
all_vars = as.data.frame(all_vars)

table(all_vars$Disease)

## Filter on readcounts (sanity check):
all_vars_RCfiltered = subset(all_vars, subset=(nalt >= 5 & (talt >= 5 | (is.na(talt)))))

## Filter on VAF (sanity check):
all_vars_RCfiltered_VAFfiltered = subset(all_vars_RCfiltered, subset=(nvaf >= 0.2 & (tvaf >= 0.2 | is.na(tvaf))))

dim(all_vars_RCfiltered_VAFfiltered)

table(all_vars_RCfiltered_VAFfiltered$Disease)

#### Get cohort frequency
all_vars_RCfiltered_VAFfiltered$mut_id = paste0(all_vars_RCfiltered_VAFfiltered$HUGO_Symbol,"_",all_vars_RCfiltered_VAFfiltered$HGVSg)
all_vars_RCfiltered_VAFfiltered$Allele_Frequency_class[ all_vars_RCfiltered_VAFfiltered$Allele_Frequency <= 0.0005 ] = "Rare (MAF =< 0.05%)"

vars_count = as.data.frame(table(all_vars_RCfiltered_VAFfiltered[,c("mut_id","HGVSc_short","HGVSp_short","Allele_Frequency_class")]))
vars_count = vars_count[vars_count$Freq > 0,]
colnames(vars_count) = c("Variant", "HGVSc", "HGVSp", "Allele_Frequency_class", "Number_of_cases")
vars_count$Variant = as.character(vars_count$Variant)
vars_count$Gene_name = unlist( lapply( strsplit( vars_count$Variant, "_" ), "[", 1))

write.table(vars_count, paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/varSummary_", batch_date, ".tsv"), quote = F, row.names = F, sep = "\t")
write.xlsx(vars_count, paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/varSummary_", batch_date, ".xlsx"), quote = F, rowNames = F)

for (i in unique(vars_count$Variant)){
  gene = unique(vars_count[vars_count$Variant == i, "Gene_name"])
  hgvsc = unique(vars_count[vars_count$Variant == i, "HGVSc"])
  hgvsp = unique(vars_count[vars_count$Variant == i, "HGVSp"])
  if (hgvsp == "p."){
    vars_count[ vars_count$Variant == i, "Variant2" ] = paste0(gene,"_",hgvsc)
  } else {
    vars_count[ vars_count$Variant == i, "Variant2" ] = paste0(gene,"_",hgvsp)
  }
}

p = ggplot(vars_count, aes(x=Variant2, y=Number_of_cases, group=Gene_name)) 
p = p + geom_bar(stat = "identity")
p = p + scale_color_brewer( palette = "Set1")
p = p + facet_grid(cols = vars(Allele_Frequency_class), scales = "free", space = "free") + theme(panel.margin = unit(0, "lines"), strip.text = element_blank())
p = p + ylab("Number of cases") + xlab("Variant")
p = p + theme_bw()
p = p + theme(text=element_text(size=8), axis.text.x = element_text(colour="black", size=8, angle = 90, vjust = 1, hjust = 1))
p = p + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), 
               axis.line = element_line (color = "black"), panel.background = element_blank(), 
               axis.title = element_text (face = "bold"),
               strip.text = element_text (face = "bold.italic"))


ggsave(filename = paste0(work_dir, "/5.CohortFrequencies/", batch_date, "/varSummary_allUpTo_", batch_date, ".pdf"), plot = p, device = "pdf", dpi = 300, w=9, h=6)

## JITTER PLOT
# get counts of every possible gene-type-class combination:
gene.type.class.counts = as.data.frame(table(all_vars_RCfiltered_VAFfiltered[,c("HUGO_Symbol", "Disease","Variant_Classification", "Overall_Classification","Allele_Frequency_class")]))

for (row in 1:nrow(gene.type.class.counts)){
  if (gene.type.class.counts[row,6] == 0){
    gene.type.class.counts[row,6] = NA
  }
}

# dot plot
gene.type.class.counts = gene.type.class.counts[!is.na(gene.type.class.counts$Freq),]

p = ggplot(gene.type.class.counts, aes(x=HUGO_Symbol, y=Variant_Classification, size=Freq, color = Overall_Classification, group=HUGO_Symbol)) 
p = p + geom_point(position=position_dodge2(width=1, preserve = "total"), aes(stroke=0.5))
p = p + scale_color_brewer( palette = "Set1")
p = p + facet_grid(cols = vars(Allele_Frequency_class), rows=vars(Disease), scales = "free", space = "free") + theme(panel.margin = unit(0, "lines"), strip.text = element_blank())
p = p + scale_size(range=c(4,6), breaks = c(1,2)) + theme_bw()
p = p + theme(text=element_text(size=15), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust = 0.2, hjust = 0.95))
p = p + xlab("Gene") + ylab("Variant Type") + guides(color=guide_legend("Variant Classification"), size=guide_legend("Affected patients"))
p = p + theme(legend.position = "bottom" )
p = p + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), 
               axis.line = element_line (color = "black"), panel.background = element_blank(), 
               axis.title = element_text (face = "bold"), strip.background = element_blank(), 
               strip.text = element_text (face = "bold.italic"))


ggsave(filename = paste0(work_dir, "/5.CohortFrequencies/", batch_date, "/varSummary_dotPlot_allUpTo_", batch_date, ".pdf"), plot = p, device = "pdf", dpi = 300, w=9, h=6)

#### Get list of variants to remove

### THIS CHANGES AS SAMPLE NUMBER INCREASES:
## Because the cohort size in this project was incrementaly increasing, but still small, I flagged variants
## if their frequency in the cohort becomes >1% (number of carriers / # total number of cases). I combine all cancer types for this calculation.
case_number = length(unique(samples$X.case_id))
case_number_cutoff = max(1, floor(case_number * 0.01))
vars_todel = as.character(vars_count[vars_count$Number_of_cases >= case_number_cutoff,"Variant"])

## Tag cohort level common variants
all_vars_RCfiltered_VAFfiltered$Cohort_frequency_notes[ all_vars_RCfiltered_VAFfiltered$mut_id %in% vars_todel ] = paste0("Variant reported in >= ", case_number_cutoff, " out of ", case_number, " PECGS cases as of ", batch_date, " (Cohort frequency >= 1%); Reporting is at the discretion of the genetic counselor.")
all_vars_RCfiltered_VAFfiltered$Cohort_frequency_notes[ !(all_vars_RCfiltered_VAFfiltered$mut_id %in% vars_todel) ] = paste0("Variant reported in < ", case_number_cutoff, " out of ", case_number, " PECGS cases as of ", batch_date, " (Cohort frequency < 1%); Reporting is at the discretion of the genetic counselor.")

## Write out per case outputs - disease specific and ACMG SF; only new cases

cols = c("Case_ID", "Disease", "HUGO_Symbol", "Chromosome", "Start", "Stop",
         "Reference", "Alternate", "Variant_Classification",
         "HGVSg", "HGVSc", "HGVSc_short","HGVSp", "HGVSp_short", 
         "Allele_Frequency", "Allele_Frequency_class", "nref", "nalt", "nvaf", "tref", "talt", "tvaf", "VEP_Most_Severe_Consequence", 
         "Positive_Evidence", "Negative_Evidence", "Positive_CharGer_Score", "Negative_CharGer_Score",
         "CharGer_Score", "ClinVar_Pathogenicity", "ACMG_Classification", "CharGer_Classification", 
         "Overall_Classification", "ClinVar_Traits","CharGer_Summary",
         "Genotype", "X1KGenomes_AF", "ExAC_AF", "gnomAD_AF", "Gene_type",
         "ACMG_SF_gene", "PanCancer_gene", "MM_specific_gene", "CHOL_specific_gene", "CRC_specific_gene", 
         "Manual_review", "Notes", "Cohort_frequency_notes")

all_vars_RCfiltered_VAFfiltered = unique(all_vars_RCfiltered_VAFfiltered)

for (case in current_cases){
  out_dir = paste0(work_dir, "/5.CohortFrequencies/", batch_date, "/outputs/", case)
  dir.create( out_dir, recursive = T )
  out = data.frame(matrix(ncol = length(cols), nrow = 0))
  colnames(out) = cols
  write.table(out, paste0(out_dir,"/",case, "_charged_reviewed_germline_ACMG_SF.tsv"), row.names = F, quote = F, sep = "\t" )
  write.table(out, paste0(out_dir,"/",case, "_charged_reviewed_germline_RARE.tsv"), row.names = F, quote = F, sep = "\t" )
  write.table(out, paste0(out_dir,"/",case, "_charged_reviewed_germline_LOW_FREQ.tsv"), row.names = F, quote = F, sep = "\t" )
}

for (case in current_cases){
  if (case %in% unique(all_vars_RCfiltered_VAFfiltered$Case_ID)){
    out = all_vars_RCfiltered_VAFfiltered[all_vars_RCfiltered_VAFfiltered$Case_ID == case,]
    cancer_type = unique(out$Disease)
    out = out[,cols]
    out_acmg = subset(out, subset = (ACMG_SF_gene == "yes")) # report acmg results as separate file
    out_dir = paste0(work_dir, "/5.CohortFrequencies/", batch_date, "/outputs/", case)
    dir.create(out_dir, recursive = T )
    if (cancer_type == "MM"){
      out_rare = subset(out, subset = (Allele_Frequency <= 0.0005 & (PanCancer_gene == "yes" | MM_specific_gene == "yes")))
      out_low_freq = subset(out, subset = (Allele_Frequency > 0.0005 & Allele_Frequency < 0.01 & (PanCancer_gene == "yes" | MM_specific_gene == "yes")))
      write.table(out_acmg, paste0(out_dir,"/",case, "_charged_reviewed_germline_ACMG_SF.tsv"), row.names = F, quote = F, sep = "\t" )
      write.table(out_rare, paste0(out_dir,"/",case, "_charged_reviewed_germline_RARE.tsv"), row.names = F, quote = F, sep = "\t" )
      write.table(out_low_freq, paste0(out_dir,"/",case, "_charged_reviewed_germline_LOW_FREQ.tsv"), row.names = F, quote = F, sep = "\t" )
    } else if (cancer_type == "CHOL"){
      out_rare = subset(out, subset = (Allele_Frequency <= 0.0005 & (PanCancer_gene == "yes" | CHOL_specific_gene == "yes")))
      out_low_freq = subset(out, subset = (Allele_Frequency > 0.0005 & Allele_Frequency < 0.01 & (PanCancer_gene == "yes" | CHOL_specific_gene == "yes")))
      write.table(out_acmg, paste0(out_dir,"/",case, "_charged_reviewed_germline_ACMG_SF.tsv"), row.names = F, quote = F, sep = "\t" )
      write.table(out_rare, paste0(out_dir,"/",case, "_charged_reviewed_germline_RARE.tsv"), row.names = F, quote = F, sep = "\t" )
      write.table(out_low_freq, paste0(out_dir,"/",case, "_charged_reviewed_germline_LOW_FREQ.tsv"), row.names = F, quote = F, sep = "\t" )
    } else if (cancer_type == "CRC"){
      out_rare = subset(out, subset = (Allele_Frequency <= 0.0005 & (PanCancer_gene == "yes" | CRC_specific_gene == "yes")))
      out_low_freq = subset(out, subset = (Allele_Frequency > 0.0005 & Allele_Frequency < 0.01 & (PanCancer_gene == "yes" | CRC_specific_gene == "yes")))
      write.table(out_acmg, paste0(out_dir,"/",case, "_charged_reviewed_germline_ACMG_SF.tsv"), row.names = F, quote = F, sep = "\t" )
      write.table(out_rare, paste0(out_dir,"/",case, "_charged_reviewed_germline_RARE.tsv"), row.names = F, quote = F, sep = "\t" )
      write.table(out_low_freq, paste0(out_dir,"/",case, "_charged_reviewed_germline_LOW_FREQ.tsv"), row.names = F, quote = F, sep = "\t" )
    } 
  }
}

### Copy new files to PECGS Results folder (for Jay)
current_cases_2cancer = samples[samples$X.case_id %in% current_cases,]

for (case in current_cases){
  print(case)
  cancer_type = unique(current_cases_2cancer[current_cases_2cancer$X.case_id == case, "cancer_type"])
  out_dir2 = paste0(work_dir,"/TEST_RESULTS/", cancer_type, "/", case, "/pipeline_reports/germline")
  print(out_dir2)
  dir.create( out_dir2, recursive = T)
  files = list.files(paste0(work_dir, "/5.CohortFrequencies/",batch_date,"/outputs"), recursive = T, pattern = case, full.names = T)
  file.copy(files, out_dir2, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
}

## Write out master table

all_vars_RCfiltered_VAFfiltered = all_vars_RCfiltered_VAFfiltered[order(all_vars_RCfiltered_VAFfiltered$Case_ID),]

## Remove variants in genes non-specific to cancer type

all_vars_RCfiltered_VAFfiltered = all_vars_RCfiltered_VAFfiltered[,cols]

all_vars_RCfiltered_VAFfiltered_mm = all_vars_RCfiltered_VAFfiltered[all_vars_RCfiltered_VAFfiltered$Disease == "MM",]
all_vars_RCfiltered_VAFfiltered_mm = subset(all_vars_RCfiltered_VAFfiltered_mm, subset = (ACMG_SF_gene == "yes" | PanCancer_gene == "yes" | MM_specific_gene == "yes"))

all_vars_RCfiltered_VAFfiltered_chol = all_vars_RCfiltered_VAFfiltered[all_vars_RCfiltered_VAFfiltered$Disease == "CHOL",]
all_vars_RCfiltered_VAFfiltered_chol = subset(all_vars_RCfiltered_VAFfiltered_chol, subset = (ACMG_SF_gene == "yes" | PanCancer_gene == "yes" | CHOL_specific_gene == "yes"))

all_vars_RCfiltered_VAFfiltered_crc = all_vars_RCfiltered_VAFfiltered[all_vars_RCfiltered_VAFfiltered$Disease == "CRC",]
all_vars_RCfiltered_VAFfiltered_crc = subset(all_vars_RCfiltered_VAFfiltered_crc, subset = (ACMG_SF_gene == "yes" | PanCancer_gene == "yes" | CRC_specific_gene == "yes"))

all_vars_RCfiltered_VAFfiltered_final = rbind(all_vars_RCfiltered_VAFfiltered_mm, all_vars_RCfiltered_VAFfiltered_chol, all_vars_RCfiltered_VAFfiltered_crc)
all_vars_RCfiltered_VAFfiltered_final = all_vars_RCfiltered_VAFfiltered_final[order(all_vars_RCfiltered_VAFfiltered_final$Case_ID),]

all_vars_RCfiltered_VAFfiltered_final$Allele_Frequency_class = gsub("≤", "<=", all_vars_RCfiltered_VAFfiltered_final$Allele_Frequency_class)
all_vars_RCfiltered_VAFfiltered_final$Cohort_frequency_notes = gsub("≤", "<=", all_vars_RCfiltered_VAFfiltered_final$Cohort_frequency_notes)

write.table(all_vars_RCfiltered_VAFfiltered_final, paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/outputs/all_reported_germline_variants_v", batch_date, ".tsv"), quote = F, row.names = F, sep = "\t")
write.xlsx(all_vars_RCfiltered_VAFfiltered_final, paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/outputs/all_reported_germline_variants_v", batch_date, ".xlsx"), quote = F, rowNames = F)

#### Create symbolic links to files
file.remove("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview_TEST/TEST_RESULTS/all_reported_germline_variants_current.tsv")
file.symlink(paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/outputs/all_reported_germline_variants_v", batch_date, ".tsv"), "/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview_TEST/TEST_RESULTS/all_reported_germline_variants_current.tsv")

file.remove("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview_TEST/TEST_RESULTS/all_reported_germline_variants_current.xlsx")
file.symlink(paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/outputs/all_reported_germline_variants_v", batch_date, ".xlsx"), "/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview_TEST/TEST_RESULTS/all_reported_germline_variants_current.xlsx")

file.remove("/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview_TEST/TEST_RESULTS/varSummary_current.xlsx")
file.symlink(paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/varSummary_", batch_date, ".xlsx"), "/storage1/fs1/dinglab/Active/Projects/PECGS/PECGS_analysis/GermlineReview_TEST/TEST_RESULTS/varSummary_current.xlsx")

#### Summary plots - current batch ####

## Variant summary
## Relabel variant type

path.lpath.relabeled = all_vars_RCfiltered_VAFfiltered[all_vars_RCfiltered_VAFfiltered$Case_ID %in% current_cases,]

for (i in 1:nrow(path.lpath.relabeled)) {
  if (path.lpath.relabeled$Variant_Classification[i] == "synonymous_variant") {
    path.lpath.relabeled[i,"Variant_Classification"]  =  "synonymous"
  } else {
    if (path.lpath.relabeled$Variant_Classification[i] == "stop_retained_variant") {
      path.lpath.relabeled[i,"Variant_Classification"]  =  "stop retained" 
    } else {
      if (path.lpath.relabeled$Variant_Classification[i] == "non_coding_transcript_exon_variant") {
        path.lpath.relabeled[i,"Variant_Classification"]  =  "non coding transcript exon"
      } else {
        if (path.lpath.relabeled$Variant_Classification[i] == "frameshift_variant" & path.lpath.relabeled[i,"Alternate"] == "-"){
          path.lpath.relabeled[i,"Variant_Classification"]  =  "frameshift deletion"
        } else {
          if (path.lpath.relabeled$Variant_Classification[i] == "frameshift_variant" & path.lpath.relabeled[i,"Reference"] == "-"){
            path.lpath.relabeled[i,"Variant_Classification"]  =  "frameshift insertion"
          } else {
            if (grepl("splice",path.lpath.relabeled$Variant_Classification[i]) == TRUE & path.lpath.relabeled[i,"Alternate"] == "-"){
              path.lpath.relabeled[i,"Variant_Classification"]  =  "splice site deletion"
            } else {
              if (grepl("splice",path.lpath.relabeled$Variant_Classification[i]) == TRUE & path.lpath.relabeled[i,"Reference"] == "-"){
                path.lpath.relabeled[i,"Variant_Classification"]  =  "splice site insertion"
              } else{
                if (grepl("splice",path.lpath.relabeled$Variant_Classification[i]) == TRUE & path.lpath.relabeled[i,"Reference"] != "-" & path.lpath.relabeled[i,9] != "Alternate"){
                  path.lpath.relabeled[i,"Variant_Classification"]  =  "splice site"
                } else {
                  if (path.lpath.relabeled$Variant_Classification[i] == "stop_gained"){
                    path.lpath.relabeled[i,"Variant_Classification"]  =  "nonsense"
                  } else {
                    if (path.lpath.relabeled$Variant_Classification[i] == "stop_lost") {
                      path.lpath.relabeled[i,"Variant_Classification"]  =  "nonstop"
                    } else {
                      if (path.lpath.relabeled$Variant_Classification[i] == "missense_variant"){
                        path.lpath.relabeled[i,"Variant_Classification"]  =  "missense"
                      } else {
                        if (path.lpath.relabeled$Variant_Classification[i] == "start_lost"){
                          path.lpath.relabeled[i,"Variant_Classification"]  =  "start lost"
                        } else {
                          if (path.lpath.relabeled$Variant_Classification[i] == "inframe_deletion"){
                            path.lpath.relabeled[i,"Variant_Classification"]  =  "inframe deletion"
                          } else {
                            if (path.lpath.relabeled$Variant_Classification[i] == "inframe_insertion"){
                              path.lpath.relabeled[i,"Variant_Classification"]  =  "inframe insertion"
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


## JITTER PLOT

# get counts of every possible gene-type-class combination:
gene.type.class.counts = as.data.frame(table(path.lpath.relabeled[,c("HUGO_Symbol", "Variant_Classification", "Overall_Classification","Allele_Frequency_class")]))

for (row in 1:nrow(gene.type.class.counts)){
  if (gene.type.class.counts[row,5] == 0){
    gene.type.class.counts[row,5] = NA
  }
}

# dot plot

gene.type.class.counts = gene.type.class.counts[!is.na(gene.type.class.counts$Freq),]

p = ggplot(gene.type.class.counts, aes(x=HUGO_Symbol, y=Variant_Classification, size=Freq, color = Overall_Classification, group=HUGO_Symbol)) 
p = p + geom_point(position=position_dodge2(width=1, preserve = "total"), aes(stroke=0.5))
p = p + scale_color_brewer( palette = "Set1")
p = p + facet_grid(cols = vars(Allele_Frequency_class), scales = "free", space = "free") + theme(panel.margin = unit(0, "lines"), strip.text = element_blank())
p = p + scale_size(range=c(4,6), breaks = c(1,2)) + theme_bw()
p = p + theme(text=element_text(size=15), axis.text.x = element_text(colour="black", size=10, angle = 90, vjust = 0.2, hjust = 0.95))
p = p + xlab("Gene") + ylab("Variant Type") + guides(color=guide_legend("Variant Classification"), size=guide_legend("Affected patients"))
p = p + theme(legend.position = "bottom" )
p = p + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), 
               axis.line = element_line (color = "black"), panel.background = element_blank(), 
               axis.title = element_text (face = "bold"), strip.background = element_blank(), 
               strip.text = element_text (face = "bold.italic"))

ggsave(filename = paste0(work_dir,"/5.CohortFrequencies/", batch_date, "/outputs/varSummary_v", batch_date, ".pdf"), plot = p, device = "pdf", dpi = 300, w=9, h=4, useDingbats = F)

sessionInfo()

