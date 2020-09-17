"""
    readvcf("vcffile.vcf.gz")

Read a vcf file and return a vcf data array of all lines.
The vcf file may be compressed (bgzip) or uncompressed.
"""
function readbigvcf(vcf)

    isfile(vcf) || error("$vcf does not exist")
    bn = splitext(vcf)
    n = bn[1]

    if occursin(r".gz$", vcf)
        f = bgzip_file(vcf)
        d = readlines(f);
        bgzip_file("$n")
    else
        d = readlines(vcf);
    end

    return d;

end
