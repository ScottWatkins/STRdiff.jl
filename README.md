# STRdiff

##STRdiff

STRdiff is a project specific tool for inferring the parent-of-origin
for de novo STR mutations found in the CEPH pedigrees.

STR mutations are first identified in the parents using a tool
like HipSTR.

The program uses beagle to infer phase of the trio containing
the parent with the de novo STR and both grandparents.  The
rest of the family is then phased independently.  The new
STR allele is identified on phased chromosomes in the
children.  Because WGS data contains some genotyping issues,
possible phasing issues, new SNP mutations, and raw variants,
a consensus chromosome carring the de novo STR allele is created.
This consensus is then compared to the grandparents to identify
the most likely grandparent transmitting the de novo allele.
Where possible, the size difference is determined based on the
information in the vcf file.


[![Build Status](https://travis-ci.com/scott.watkins@genetics.utah.edu/STRdiff.jl.svg?branch=master)](https://travis-ci.com/scott.watkins@genetics.utah.edu/STRdiff.jl)
