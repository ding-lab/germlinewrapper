# germlinewrapper
detect germline variants from normal samples

Song Cao

GermineWrapper pipeline is a fully automated and modular software package designed for detection of germline variants from normal exome data. It works on LSF job scheduler. Multiple standard variant callings are included in the pipeline such as varscan, gatk and pindel.

perl germline_calling_v1.1.pl  --srg --step --sre --rdir --ref --log --q

rdir = full path of the folder holding files for this sequence run (user must provide)

srg = bam having read group or not: 1, yes and 0, no (default 1)

log = full path of the folder for saving log file; usually upper folder of rdir 

sre = re-run: 1, yes and 0, no  (default 0)

step = run this pipeline step by step. (user must provide)

ref = the human reference: 

q = which queue for submitting job; research-hpc, ding-lab, long (default)
