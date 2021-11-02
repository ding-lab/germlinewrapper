# germlinewrapper V1.1, compute1

for HG38 reference

Detect germline variants from normal sample

### Song Cao, email: scao@wustl.edu ###

### ********you must enter the directory with germlinewrapper pipeline to submit the jobs******* ###

GermineWrapper pipeline is a fully automated and modular software package designed for detection of germline variants from normal exome data. It works on LSF job scheduler. Multiple standard variant callings are included in the pipeline such as varscan, gatk and pindel.


Usage: perl $0  --srg --step --sre --rdir --ref --log --groupname --users --q

rdir = full path of the folder holding files for this sequence run (user must provide)
srg = bam having read group or not: 1, yes and 0, no (default 1)
groupname = job group name
users = user name for job group
log = full path of the folder for saving log file; usually upper folder of rdir 
sre = re-run: 1, yes and 0, no  (default 0)
step run this pipeline step by step. (user must provide)
ref the human reference: 
q which queue for submitting job; research-hpc, ding-lab, long (default)

GDC HG38: /storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa 

run_folder = full path of the folder holding files for this sequence run
step_number run this pipeline step by step. (running the whole pipeline if step number is 0)

[1]  Run gatk

[2]  Run varscan

[3]  Run pindel

[4]  Parse pindel

[5]  filter vcf

[6]  Merge calls

[7]  VCF2MAF

[8]  Generate final maf

### see work_log_test for how to run jobs in compute1 

### Details about the implementation ###


SNV Variant Calls:

* GATK4

* VarScan 2.3.8

* SNV calls are union calls from VarScan and GATK.


Indel Variant Calls:

* Pindel

* VarScan 2.3.8

* GATK4

* Indel calls are from variants called by both GATK and VarScan or Pindel.

* Merging SNV and Indel calls

* Annotation resulting in annotated MAF

* filter large indels longer than 100 bps
