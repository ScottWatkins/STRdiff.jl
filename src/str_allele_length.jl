"""
    str_allele_length("myvcf.gz", chrpos)

Build a simple lookup table mapping the allele
to the length of that allele for a give STR or
indel variant in a vcf file.  Vcf must be compressed
and indexed; location should be given as 1:123456.

Dependencies: bcftools
"""
function str_allele_length(vcf, chrpos)

    lup = Dict{Int64, Int64}()

    if Sys.which("bcftools") == nothing
        error("bcftools not found. Please add bcftools to your path.")
    end

    run(`bcftools index -f  $vcf`)

    (chr, pos) = split(chrpos, ':')

    s = read(`bcftools view -H -r $chrpos $vcf`, String)

    if ( length(s) < 1 ) && ( occursin(r"_fam_", vcf) )
        error("the STR specified at $chrpos was not found!")
    end

    t = split(s, '\t')

    lup[0] = length(t[4]) #ref allele
    lup[-1] = 0           #missing data (./. recoded to -1/-1)

    x = split(t[5], ',')

    for (i, j) in enumerate(x)  #all alt alleles
        k = length(j)
        lup[i] = k
    end

    return lup

end
