# GermlineWrapper v2.4

HG38 germline variant calling pipeline for normal exome sequencing data  
Designed for LSF-based HPC environments (compute1)

---

## Overview

GermlineWrapper is a fully automated and modular pipeline for detecting germline variants from normal exome sequencing data. It integrates multiple widely used variant calling tools and provides a complete workflow from raw BAM files to annotated variants, clinical classification, and manual review outputs.

Supported reference: hg38 (GDC GRCh38)

---

## Author

Song Cao  
Washington University School of Medicine  
Email: scao@wustl.edu

---

## Key Features

- Multi-caller germline variant detection (GATK4, VarScan, Pindel)
- VEP annotation
- CharGer classification
- Readcount integration
- Automated report generation
- IGV manual review support
- LSF job scheduler compatible
- Modular step-by-step execution

---

## Updates in v2.4

- Implemented Fernanda VEP annotation workflow  
- Tumor readcount defaults to 0 if unavailable  
- Added downstream steps 15–18 (parse CharGer, reporting, IGV review)

---

## Usage

perl germlinewrapper.pl --rdir <run_dir> --log <log_dir> --groupname <group> --users <username> --ref <reference.fa> --step <step_number> --q <queue>

---

## Pipeline Steps

1. Run GATK  
2. Run VarScan  
3. Run Pindel  
4. Parse Pindel  
5. Filter VCF  
6. Merge calls  
7. VCF2MAF  
8. Generate final MAF  
9. Run bam-readcount  
10. Add readcount  
11. Generate MAF with readcount  
12. Generate VEP input  
13. Run VEP  
14. Run CharGer  
15. Postprocess CharGer  
16. Add readcount and IGV session  
17. Generate final report  
18. Generate IGV scripts  

---

## Contact

Song Cao  
scao@wustl.edu
