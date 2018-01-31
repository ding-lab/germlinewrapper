# germlinewrapper
detect germline variants from normal samples

Song Cao

GermineWrapper pipeline is a fully automated and modular software package designed for detection of germline variants from normal exome data. It works on LSF job scheduler. Multiple standard variant callings are included in the pipeline such as varscan, gatk and pindel.

Additional notes

* All paths with respect to host 

TODO:
set up centromeres, etc in `/image.setup/C_Centromeres/pindel-centromere-exclude.bed`
set up reference (or just use SomaticWrapper's)
confirm that picard is needed, and install in SomaticWrapper docker image
