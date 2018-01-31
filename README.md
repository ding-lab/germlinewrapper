# Germline Wrapper d2

Detect germline variants from normal samples

Song Cao and Matt Wyczalkowski

GermineWrapper pipeline is a fully automated and modular software package
designed for detection of germline variants from normal exome data. 
Designed to run in MGI and standard docker envrionments. 

Variant callers used: varscan, gatk and pindel.

## Additional notes

* All paths with respect to host 

TODO:
set up centromeres, etc in `/image.setup/C_Centromeres/pindel-centromere-exclude.bed`
set up reference (or just use SomaticWrapper's)
confirm that picard is needed, and install in SomaticWrapper docker image

* GermlineWrapper will be installed into SomaticWrapper docker image
* Add picard to SW docker image

