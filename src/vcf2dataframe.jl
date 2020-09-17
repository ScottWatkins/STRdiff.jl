"""
    vcf2dataframe(vcf_data_array)

Convert a vcf data array from the readbigvcf.jl function
and return a standard data frame.

Requires: DataFrames
"""
function vcf2dataframe(vcf_data_array)

vcf = vcf_data_array
hapmat = []
header = []
dataline = []
lup = Dict()

for i in eachindex(vcf)

    if occursin(r"^##", vcf[i])

        continue

    elseif occursin(r"^#CHROM", vcf[i])

        ids = split(vcf[i], r"\t")
        ids[1] = replace(ids[1], "#" => "")

        push!(header, ids[1], ids[2])

        for a = 10:length(ids)
            push!(header, ids[a] * "p", ids[a] * "q")
        end

        hapmat = [hapmat; header]

    else

        g = split(vcf[i], '\t')

        push!(dataline, g[1], g[2])

        for i = 10:length(g)
            a = split(g[i], r"[\||\:|\/]")
            push!(dataline, a[1], a[2])
        end

        dataline = replace.(dataline, "." => "-1") #missing to -1

        addbases = false #routine to add bases calls

        if addbases == true
            alts = enumerate(split(g[5], ","))

            if (length(g[4]) +  length(g[5])) > 2 #indels remain numeric
                continue
            end

            dataline = replace.(dataline, r"0" => g[4])

            for (i,v) = alts
                dataline = replace.(dataline, Regex("$i") => v)
            end
        end

        hapmat = [hapmat dataline]

        dataline = []
    end
end

p = permutedims(hapmat)
pm = p[2:end, :]
h = p[1, :]
t = DataFrame(pm)
tt = rename!(t, Symbol.(h))

return tt

end
