# germlinewrapper V1.1

for HG38 reference

Detect germline variants from normal sample

### Song Cao, email: scao@wustl.edu ###

### ********you must enter the directory with germlinewrapper pipeline to submit the jobs******* ###

GermineWrapper pipeline is a fully automated and modular software package designed for detection of germline variants from normal exome data. It works on LSF job scheduler. Multiple standard variant callings are included in the pipeline such as varscan, gatk and pindel.

perl germlinewrapper.pl  --srg --step --sre --rdir --ref --log --q

rdir = full path of the folder holding files for this sequence run (user must provide)

srg = bam having read group or not: 1, yes and 0, no (default 1)

log = full path of the folder for saving log file; usually upper folder of rdir 

sre = re-run: 1, yes and 0, no  (default 0)

step = run this pipeline step by step. (user must provide)

ref = the human reference: 

q = which queue for submitting job; research-hpc, ding-lab, long (default)

### An example for running step on MGI cluster ###

perl germlinewrapper.pl --rdir /gscmnt/gc2521/dinglab/scao/cptac3/hg38/germline_unfinished --ref /gscmnt/gc2521/dinglab/mwyczalk/somatic-wrapper-data/image.data/A_Reference/GRCh38.d1.vd1.fa --q research-hpc --log /gscmnt/gc2521/dinglab/scao/cptac3/hg38 --step 1

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
