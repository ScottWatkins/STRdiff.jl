# STRdiff

STRdiff is a project specific tool for inferring the parent-of-origin for de novo STR mutations found in the CEPH pedigrees.

STR mutations are first identified in the parents using a tool like HipSTR.

The program uses beagle to infer chromosome phase in the trio containing the parent with the de novo STR and both grandparents. The rest of the family is then also phased independently. The new STR allele is then identified on phased chromosomes in the children. Because phased WGS data will contain some genotyping issues, possible phasing issues, new SNP mutations, and rare variants, a consensus chromosome carrying the de novo STR allele is created. This consensus is then compared to the grandparents to identify the most likely grandparental chromosome transmitting the de novo allele. Where possible, the size difference is determined based on the information in the vcf file.


[![Build Status](https://travis-ci.com/scott.watkins@genetics.utah.edu/STRdiff.jl.svg?branch=master)](https://travis-ci.com/scott.watkins@genetics.utah.edu/STRdiff.jl)
