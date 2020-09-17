
include("bgzip_file.jl")
include("readbigvcf.jl")
include("vcf2dataframe.jl")
include("get_pedfile_trios.jl")
include("str_allele_length.jl")

"""
    (t, dt, p, nd, mutCon) = STR_iht(vcffile.vcf.gz, pedigree, STR_location, iid)

Read a small **phased** vcf data array and return dataframes containing 1)
transmitting chromosomes and the STR allele differences and 2) a dataframe
of all phased individuals.

Individuals in the pedigree must match those in the vcf file.  The pedigree
file may contain multiple pedigrees.

STR location format is 1:258132133

Requires: DataFrames, CSV
"""
function STR_iht(vcf, pedfile, str, iid, mutSTR, othpar)

    vcf_data_array = readbigvcf(vcf)

    hap_df = vcf2dataframe(vcf_data_array)

    (trios, peds) = get_pedfile_trios(pedfile)

    ftdf = DataFrame(ped=Any[], chr=Any[], pos=Any[], parent_chr=Any[],
    parent_sex=Any[], parent_str=Any[], parent_bp=Any[], child_chr=Any[],
    child_str=Any[], child_bp=Any[], str_diff=Any[], str_bp_diff=Any[] )

    fftdf = DataFrame(ped=Any[], chr=Any[], pos=Any[], parent_chr=Any[],
    parent_sex=Any[], parent_str=Any[], parent_bp=Any[], child_chr=Any[],
    child_str=Any[], child_bp=Any[], str_diff=Any[], str_bp_diff=Any[] )

    strchr, strpos = split(str, ':')
    stridx = argmax((hap_df[:, :POS] .== strpos))    #idx of str site

    str_sizes = str_allele_length(vcf, str)


for i in peds

    pid = i; header = []; tp = []

    t_all = get(trios, i, 0)

    for i in t_all
        line = "-------------------------\n"

        ids = split(i, ',')

#print(line, ids, line, "\n")

        for a in eachindex(ids)
            push!(header, Symbol(ids[a] * "p"), Symbol(ids[a] * "q"))
        end

#println(header)

        tdf = hap_df[ :, header]  #df of each child,father,mother trio

        c1f1 = sum(tdf[!,1] .== tdf[!,3]) / size(tdf, 1) #fraction matching
        c1f2 = sum(tdf[!,1] .== tdf[!,4]) / size(tdf, 1) #across all sites
        c1m1 = sum(tdf[!,1] .== tdf[!,5]) / size(tdf, 1)
        c1m2 = sum(tdf[!,1] .== tdf[!,6]) / size(tdf, 1)
        c2f1 = sum(tdf[!,2] .== tdf[!,3]) / size(tdf, 1)
        c2f2 = sum(tdf[!,2] .== tdf[!,4]) / size(tdf, 1)
        c2m1 = sum(tdf[!,2] .== tdf[!,5]) / size(tdf, 1)
        c2m2 = sum(tdf[!,2] .== tdf[!,6]) / size(tdf, 1)

        c1altidx = findall(tdf[!,1] .!= "0")  #alt alleles matching
        c2altidx = findall(tdf[!,2] .!= "0")  #across sites, >1 allele
        f1altidx = findall(tdf[!,3] .!= "0")  #allowed
        f2altidx = findall(tdf[!,4] .!= "0")
        m1altidx = findall(tdf[!,5] .!= "0")
        m2altidx = findall(tdf[!,6] .!= "0")

        c1f1alt = length(intersect( Set(c1altidx), Set(f1altidx) ) ) #alt allele
        c1f2alt = length(intersect( Set(c1altidx), Set(f2altidx) ) ) #matched
        c1m1alt = length(intersect( Set(c1altidx), Set(m1altidx) ) ) #counts
        c1m2alt = length(intersect( Set(c1altidx), Set(m2altidx) ) )
        c2f1alt = length(intersect( Set(c2altidx), Set(f1altidx) ) )
        c2f2alt = length(intersect( Set(c2altidx), Set(f2altidx) ) )
        c2m1alt = length(intersect( Set(c2altidx), Set(m1altidx) ) )
        c2m2alt = length(intersect( Set(c2altidx), Set(m2altidx) ) )

        c1pfm =  [c1f1, c1f2, c1m1, c1m2]  #fractions for data above
        c2pfm =  [c2f1, c2f2, c2m1, c2m2]

        c1altmatch = [c1f1alt, c1f2alt, c1m1alt, c1m2alt] ./ size(tdf, 1)
        c2altmatch = [c2f1alt, c2f2alt, c2m1alt, c2m2alt] ./ size(tdf, 1)

#println("c1pfm=", c1pfm)
#println("c1altmatch=", c1altmatch)
#println("c2pfm=", c2pfm)
#println("c2altmatch=", c2altmatch)

        chriht = ["f1", "f2", "m1", "m2"]

        #pick all max indices for 1) total fraction of sites matched
        #and total alt alleles matched
        c1maxidx = [i for (i, x) in enumerate(c1pfm) if x == maximum(c1pfm)]
        c2maxidx = [i for (i, x) in enumerate(c2pfm) if x == maximum(c2pfm)]
        c1altmaxidx = [i for (i, x) in enumerate(c1altmatch) if x == maximum(c1altmatch)]
        c2altmaxidx = [i for (i, x) in enumerate(c2altmatch) if x == maximum(c2altmatch)]

        c1max = chriht[c1maxidx]   #translated idx to f1, f2, ...
        c2max = chriht[c2maxidx]

        c1altmax = chriht[c1altmaxidx]
        c2altmax = chriht[c2altmaxidx]

#println("c1altmax=", c1altmax)
#println("c2altmax=", c2altmax)


        if maximum(c1pfm) < 0.80     #declare phasing error for >20% mismatch
            c1max = ["pe"]; c1altmax = ["pe"]
        end

        if maximum(c2pfm) < 0.80
            c2max = ["pe"]; c2altmax = ["pe"]
        end

#xxx = maximum(c1pfm)
#yyy = maximum(c2pfm)
#println("argmax_c1= ", xxx)
#println("argmax_c2= ", yyy)
#println("c1max=", c1max)
#println("c2max=", c2max)

        c1_i = intersect(c1max, c1altmax)  #take isect of best frac & best alts
        c2_i = intersect(c2max, c2altmax)  #->most likely inheritance

#println("i1max=", c1_i)
#println("i2max=", c2_i)

        function pickone(c)  #if f1,f2 or m1,m2 equal probs then choose randomly
            if length(c) == 2 && ( c[1][1] == c[2][1] )
                nc = rand(c, 1)
                return nc
            else
                return c
            end
        end

        c1_i = pickone(c1_i)
        c2_i = pickone(c2_i)


        if length(c1_i) == 1  #logic to test other parent if one parent works
            if sum(occursin.(r"pe", c1_i)) > 0
                c1_i = ["nd"]
            elseif sum(occursin.(r"f", c1_i)) > 0
                c2_i = c2_i[occursin.(r"m", c2_i), :]
            elseif sum(occursin.(r"m", c1_i)) > 0
                c2_i = c2_i[occursin.(r"f", c2_i), :]
            end
        end

#println(c1_i)
#println(c2_i)

        if length(c2_i) == 1
            if sum(occursin.(r"pe", c2_i)) > 0
                c2_i = ["nd"]
            elseif sum(occursin.(r"f", c2_i)) > 0
                c1_i = c1_i[occursin.(r"m", c1_i), :]
            elseif sum(occursin.(r"m", c2_i)) > 0
                c1_i = c1_i[occursin.(r"f", c1_i), :]
            end
        end

        if ( (length(c1_i) == 0) || (length(c1_i) ≥ 2) )
            c1_i = []
            push!(c1_i, "nd")
        end

        if ( (length(c2_i) == 0) || (length(c2_i) ≥ 2) )
            c2_i = []
            push!(c2_i, "nd")
        end

#println(c1_i)
#println(c2_i)

        h = String.(header)

        function mapid(c, h)  #translate back to ids
            k = replace.(c, "f1" => h[3])
            l = replace.(k, "f2" => h[4])
            m = replace.(l, "m1" => h[5])
            n = replace.(m, "m2" => h[6])
            return n
        end

        c1_bestids = mapid(c1_i, h)
        c2_bestids = mapid(c2_i, h)

#println(c1_bestids)
#println(c2_bestids)

        #flatten arrays and push results to a df
        c1_i = join(reshape(c1_i, (1, length(c1_i) ) ), ",")
        c2_i = join(reshape(c2_i, (1, length(c2_i) ) ), ",")
        c1_bestids = join(reshape(c1_bestids, (1, length(c1_bestids) ) ), ",")
        c2_bestids = join(reshape(c2_bestids, (1, length(c2_bestids) ) ), ",")


        print_denovo_trio = true
        if print_denovo_trio == true
            if occursin(Regex("$iid[p|q]"), names(tdf)[1]) #print only denovo parent trio
                println("INFO: phased chromosomes for trio with the de novo STR.")
                println("INFO: de novo STR site genotypes are -1 (unphasable)")
                println(tdf)
                println("$iid p transmitted from: ", c1_bestids,
                "\n", "$iid q transmitted from: ", c2_bestids, "\n")
            end
        end

        c1chr = h[1]
        c2chr = h[2]

        if (c1_bestids != "nd")
            p1_allele = parse(Int64, (tdf[ stridx, "$c1_bestids" ] ) )
            c1_allele = parse(Int64, (tdf[ stridx, "$c1chr" ] ) )
            p1_length = str_sizes[p1_allele]
            c1_length = str_sizes[c1_allele]
            c1strdiff = abs(p1_allele - c1_allele)
            c1bpdiff  = p1_length - c1_length
        else
            c1strdiff = "nd"; c1_allele = "nd"; c1_length = "nd"
            c1bpdiff  = "nd"; p1_allele = "nd"; p1_length = "nd"
        end

        if (c2_bestids != "nd")
            p2_allele = parse(Int64, (tdf[ stridx, "$c2_bestids" ] ) )
            c2_allele = parse(Int64, (tdf[ stridx, "$c2chr" ] ) )
            p2_length = str_sizes[p2_allele]
            c2_length = str_sizes[c2_allele]
            c2strdiff = abs(p2_allele - c2_allele)
            c2bpdiff  = p2_length - c2_length
        else
            c2strdiff = "nd"; c2_allele = "nd"; c2_length = "nd"
            c2bpdiff  = "nd"; p2_allele = "nd"; p2_length = "nd"
        end

        if occursin(r"f", c1_i)
            p1_sex = "male"
        elseif occursin(r"m", c1_i)
            p1_sex = "female"
        else
            p1_sex = "nd"
        end

        if occursin(r"f", c2_i)
            p2_sex = "male"
        elseif occursin(r"m", c2_i)
            p2_sex = "female"
        else
            p2_sex = "nd"
        end


        push!(ftdf, [pid strchr strpos c1_bestids p1_sex p1_allele p1_length c1chr c1_allele c1_length c1strdiff c1bpdiff] )
        push!(ftdf, [pid strchr strpos c2_bestids p2_sex p2_allele p2_length c2chr c2_allele c2_length c2strdiff c2bpdiff] )

        r1 = []   #reinit for next trio
        r2 = []
        header = []
    end

    fracnd = round(sum(ftdf.str_diff .== "nd") / size(ftdf, 1), digits=2) #frac bad

    fftdf = ftdf[ .&(ftdf.str_diff .!= 0, ftdf.str_diff .!= "nd"), :]

#Next, extract the rest of the family that carries the mutant STR allele
#and get the consensus haplotype with the mutation and compare
#to the grandparents later.

    mutCon = nothing
    num_mut_iht = 0

#println(ftdf)

    if size(hap_df,2) > 8  #run fam but not trio here

        mut_offspr = (ftdf[ ftdf.child_str .== mutSTR, [:child_chr]]) #df
        mut_ids = vec(Array(mut_offspr))  #array, 2209p ...
        mut_df = hap_df[:, Symbol.(mut_ids)]

#println(mut_df)
#println(mut_offspr)
#println(mut_ids)

        if sum(occursin.(Regex(othpar), mut_ids)) > 0 #other parent has de novo
            println("WARN: one or more phased chromosomes from parent $othpar
also carries the mutant STR allele.")
            println("INFO: omitting phased chromosomes containing the de novo STR
allele that match the other parent $othpar at a rate of ≥ 0.85")
            othchrid = mut_offspr[ occursin.(Regex(othpar), mut_ids), :]
            othchrid = othchrid[1,1]
            othchr = Array(mut_df[:, Symbol(othchrid)])

            p_oth = Float64[]
            not_othchr = Int64[]

            println(mut_df)

            for i in 1:size(mut_df, 2)
                p = sum(othchr .== mut_df[:,i]) / size(mut_df, 1)
                push!(p_oth, p)
            end

            for i in 1:length(p_oth)
                if p_oth[i] ≥ 0.85
                    continue
                else
                    push!(not_othchr, i)
                end
            end

            mut_df = mut_df[:, not_othchr]

        end

        if size(mut_df, 2) == 0

            println("ERROR: consensus phased chromosome with de novo STR also matches
the other parent (>90%). This result may result can be caused by 1) very high
similarity among grandparent haplotypes, 2) incorrect phasing, or 3) a false positive
de novo STR.  Try rerunning this site with different haplotype lengths\n")

            num_mut_iht = -9
            return(fftdf, ftdf, hap_df, fracnd, mutCon, str_sizes, num_mut_iht)
        end

        av = sum(eachcol(mut_df .== "1")) / size(mut_df,2) #generate a consensus
        mutCon = Array{String,1}(undef, length(av) )       #haplotype (mutCon)

        for i in eachindex(av)
            if av[i] ≥ 0.5
                mutCon[i] = "1"
            else
                mutCon[i] = "0"
            end
        end

        num_mut_iht = size(mut_df,2)

        println("INFO: Found $num_mut_iht offspring carrying the de novo allele $mutSTR\n")
        println("TABLE: Phased offspring chromosomes with de novo allele:")
        println(mut_df)

    end

    sort!(fftdf, :parent_chr)

    return (fftdf, ftdf, hap_df, fracnd, mutCon, str_sizes, num_mut_iht)

end
end
