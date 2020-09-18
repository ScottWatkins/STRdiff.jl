"""
    getSTRdiff(datafile, haplotype_size;
        show="full_table|full_pedigrees", showhaps=false)

STRdiff is a program to extract and phase the region surrounding a suspected
de novo STR event in the CEPH pedigrees.

All trios from a pedigree are processed to identify any difference in the
STR allele length between grandparent and parent.

A table of these differences, based on local chromosomal phasing, and the sex
of the transmitting parent is returned.

This version requires 1) that the pedigree has two parents and grandparents
and 2) that the de novo event occurred in one of the CEPH parents.

Inheritance of the new STR allele found on the phased chromosomes of the
CEPH children is compared with the phased chromosomes in the grandparents to
determine which grandparental chromosome transmitted the mutation.

It is strongly recommended that a given STR be analyzed at different
haplotype lengths and that only runs with a unique solution score of 
≥5 be considered for additional analysis.  Use the repSTRdiff function
to automate multiple runs 

Dependencies: bcftools, bgzip

Input file format example: (single spaced, tab-delimited, header optional)\n
#PEDIGREE       VCF_FILE                IID     STRLOC         STR_TYPE(bp)\n
1463.ped        1463_phase_in.vcf.gz    2188    2:174994186    2\n
1337.ped        1337.vcf.gz             1295    1:2581321      4\n
"""
function getSTRdiff(infile, hapsize::Int64; kwargs...)

kw = Dict(kwargs)
sh = get(kw, :show, nothing)
s_haps = get(kw, :showhaps, false)

fname = splitext(infile)
fn = fname[1]
loc_list = String[]

if (hapsize < 1000) || (hapsize > 500000)
    error("hapsizes, bp on each side of the STR, should be between 10000 and 320000.")
end

OUT = open("$fn.$hapsize.out", "w")
header =  "#ped\tiid\tloc\ttranmitting_grandparent\ttransmitting_grandparent_sex\ttransmitting_grandparent_str\ttransmitting_grandparent_str_size\tparent_str\tparent_str_size\tstr_allele_change\tsize_change_bp\toffspring_with_de_novo_allele\tgp_match_probability\tuniq_soln_score\thapsize\tstr_type"
println(OUT, header)
close(OUT)

open(infile) do file

    for line in eachline(file)

        if occursin(r"^#", line)
            continue
        end

        l = '-'^80

        println("\n$l")
        printstyled("START RUN\nINFO: Processing $line\nINFO: Haplotype size = $hapsize\n", color=:cyan, bold=true)
        println("$l")

        z = preprocessDenovoSTRs(line, window=hapsize)

        datin = split(line, '\t')
        if length(datin) != 5
            error("expecting 5 fields for each STR location,\nfound line: $datin")
        end
        pedname = splitext(datin[1])
        ped = pedname[1]
        shf = datin[1]
        did = datin[3]
        expbp = datin[5]

        sex = hashfile("$shf", 2, [5])
        fathers = hashfile("$shf", 2, [3])

        f1 = findmax(countmap(collect(values(fathers))))[2][1]
        mothers = hashfile("$shf", 2, [4])
        m1 = findmax(countmap(collect(values(mothers))))[2][1]
        othpar = ""

        if did == f1
            othpar = m1
        elseif did == m1
            othpar = f1
        end

        println("INFO: Parent with de novo STR is $did")
        println("INFO: Other parent id is $othpar")

        fa = z[1][1]; fb = z[1][2]; site = z[1][3]; iid = z[1][4]
        ta = z[1][5]; tb = z[1][6]; mutSTR = z[1][7]; par_genos = z[1][8]

        (p, q, r, u, mutCon, fa_strs, num_iht_mut) = STR_iht(fa, fb, site, iid, mutSTR, othpar; showhaps=s_haps)
        (tp, tq, tr, tu, mutCon_empty, t_strs, num_iht_mut_null) = STR_iht(ta, tb, site, iid, mutSTR, othpar; showhaps=s_haps)

        if (mutCon != nothing) && ( size(tr, 1) != length(mutCon) )
            error("Variant site $site has failed for hapsize $hapsize.
Try a different haplotype size or consider removing this site
if it fails at multiple hapsizes.")
            exit()
        end

        if num_iht_mut == -9

            na = "na"; rr = -9
            OUT = open("$fn.$hapsize.out", "a+")
            println(OUT, "$ped\t$iid\t$site\t$na\t$na\t$na\t$na\t$na\t$na\t$na\t$na\t$num_iht_mut\t$na\t$rr\t$hapsize\t$expbp")
            close(OUT)
            println("\n$l\n$l")
            continue

        else

            match_p = []

            for i in 3:6
                p = sum(tr[!,i] .== mutCon)/length(mutCon)
                push!(match_p, p)
            end

            parent_idx = rand(findall( match_p .== maximum(match_p) ) ) + 2

            match_p_sorted = sort(match_p, rev=true)

            mdiff = round(match_p_sorted[1] - match_p_sorted[2], digits=2) * 100

            uss = Int64(floor(mdiff))

            rflag = ""

            if length(findall( match_p .== maximum(match_p) ) ) > 1
                rflag = "rflag"
                uss = 0
            end


            trans_p = names(tr,parent_idx)  #transmitting parent id

            trans_ps = trans_p[1]  #convert to plain string for patmatch

            #Test transmitting parent against parents genotypes and get
            #the alleles.  Test homs then hets

            psize = []
            pg = []

            if occursin(par_genos[1], trans_ps) && (par_genos[2] == par_genos[3])

                push!(pg, parse(Int, par_genos[2]))
                push!(psize, get(t_strs, pg[1], 0))

            elseif occursin(par_genos[4], trans_ps) && (par_genos[5] == par_genos[6])

                push!(pg, parse(Int, par_genos[5]))
                push!(psize, get(t_strs, pg[1], 0))

            elseif occursin(par_genos[1], trans_ps) && (par_genos[2] != par_genos[3])

                pg_1 = parse(Int, par_genos[2])
                pg_2 = parse(Int, par_genos[3])
                push!(pg, pg_1, pg_2)
                push!(psize, get(t_strs, pg_1, 0))
                push!(psize, get(t_strs, pg_2, 0))

            elseif occursin(par_genos[4], trans_ps) && (par_genos[5] != par_genos[6])

                pg_1 = parse(Int, par_genos[5])
                pg_2 = parse(Int, par_genos[6])
                push!(pg, pg_1, pg_2)

                push!(psize, get(t_strs, pg_1, 0))
                push!(psize, get(t_strs, pg_2, 0))

            else

                println("ERROR: unable to uniquely determine all phased STR
alleles transmitting the de novo event, please vary the hapsize value.")
                rr = -9
                na = "na"
                println("$par_genos $trans_ps")
                OUT = open("$fn.$hapsize.out", "a+")
                println(OUT, "$ped\t$iid\t$site\t$na\t$na\t$na\t$na\t$na\t$na\t$na\t$na\t$na\t$na\t$rr\t$hapsize\t$expbp")
                close(OUT)
                continue

            end

            dsize = get(fa_strs, mutSTR, 0)

            if length(psize) > 1

                dsize_diff1 = dsize - psize[1]
                dsize_diff2 = dsize - psize[2]
                str_s_diff1 = mutSTR - pg[1]
                str_s_diff2 = mutSTR - pg[2]

                if abs(dsize_diff1) == abs(dsize_diff2)
                    dsize_diff = join(["⨦", abs(dsize_diff1)], "")
                elseif abs(dsize_diff1) != abs(dsize_diff2)
                    dsize_diff = NaN
                end

                if abs(str_s_diff1) == abs(str_s_diff2)
                    str_s_diff = join(["⨦", abs(str_s_diff1)], "")
                elseif abs(str_s_diff1) != abs(str_s_diff2)
                    str_s_diff = NaN
                end

            else

                dsize_diff = dsize - psize[1]
                str_s_diff = mutSTR - pg[1]

            end

            trans_pid = trans_ps[1:end-1]
            trans_ps_sex = get(sex, trans_pid, 0)
            tps = trans_ps_sex[1]

            function namesex(x)
                tps = ""
                if x == "2"
                    tps = "female"
                elseif x == "1"
                    tps = "male"
                else
                    tps = "unknown"
                end
                return tps
            end

            if sh == "full_tables"
                println(q, "\n", tq)
            elseif sh == "full_pedigrees"
                println(r, "\n", tr)
            end

            mprobs = round.(match_p, digits=6)
            tps_sex = namesex(tps)
            println()
            println("Grandparents match probability with offspring consensus: ", round.(match_p, digits=6) )
            println("Grandparent transmitting the de novo allele: $trans_ps  $rflag")
            println("Transmitting grandparent's sex: $tps_sex")
            println("Unique solution score: $uss")
            println("Number of de novo alleles found in phased offspring chromosomes: $num_iht_mut")

            ppg = join(pg, ",")
            ppsize = join(psize, ",")

            println("Grandparental STR allele and size: $ppg and $ppsize bp")
            println("Parental de novo STR allele and size: $mutSTR and $dsize bp")
            println("De novo STR allele size difference (bp): $dsize_diff bp")

            OUT = open("$fn.$hapsize.out", "a+")
            println(OUT, "$ped\t$iid\t$site\t$trans_ps\t$tps_sex\t$ppg\t$ppsize\t$mutSTR\t$dsize\t$str_s_diff\t$dsize_diff\t$num_iht_mut\t$mprobs\t$uss\t$hapsize\t$expbp")
            close(OUT)

            println("$l\n$l")
        end

        run(`bash -c "rm -f *fam*"`)
        run(`bash -c "rm -f *trio*"`)

    end

end

end
