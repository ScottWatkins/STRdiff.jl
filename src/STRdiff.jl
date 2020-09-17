module STRdiff

#check for external dependencies: 1) bcftools, 2) bgzip
bcf = Sys.which("bcftools")
bgz = Sys.which("bgzip")

if bcf == nothing || bgz == nothing
    error("could not find bcftools or bgzip. Please make sure bcftools (ver. 1.6+) and bgzip are in your path.")
end


using CSV
using DataFrames
using DelimitedFiles
using StatsBase

export check_vcf_sites
export getSTRdiff
export repSTRdiff

include("check_vcf_sites.jl")
include("split_trio_from_ped.jl")
include("preprocessDenovoSTRs.jl")
include("STR_iht.jl")
include("hashfile.jl")
include("getSTRdiff.jl")
include("repSTRdiff.jl")

end
