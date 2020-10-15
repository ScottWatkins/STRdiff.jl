"""
    repSTRdiff("runparams.txt", [hapsize_1, hapsize_2, ... hapsize_N];
        keepall=false, showhaps=false, min_uss_score=10)

Run getSTRdiff for multiple user-specified haplotype sizes.
Retrieve the best estimates from all runs and output a final table.

See getSTRdiff for additional details on input format.

Options:

    keepall        Keep all intermediate files [false]

    showhaps       Show the offspring and de novo trio haplotypes [false]

    min_uss_score  Minimum value for the unique solution score.  For values
                   below this number, both chromosomes from each of the
                   grandparents will be tested, and the minimum difference among
                   probabilities for both male and both female grandparental
                   chromosomes will be shown parenthetically next to the uss
                   score (e.g. 3,(17)).  This allows the sex of the transmitter
                   to be determined in (rare) cases where there is high similarity 
                   between the transmitting parent's chromosomes and no one single
                   chromosome solution produced a high uss score.  For example,
                   [0.91, 0.90, 0.40, 0.73] would return (17), indicating support for
                   de novo transmission by the first grandparent (male)
                   even though the single chromosome solution is low.
                   This value is shown only when the uss score for the best estimate
                   falls below the min_uss_score setting [default = 10].

----------------------------------------------------------
Example sequence of commands for automation:
----------------------------------------------------------

1.  Run check\\_vcf\\_sites to clean input sites.

    check\\_vcf\\_sites("infile.vcf.gz", "datafile.txt"; c=4)

2.  Run repSTRdiff at a single haplotype length:

    repSTRdiff("runparams.txt", [20000])

    Remove any STR sites that fail from the input file. 

3.  Run repSTRdiff over a range of haplotype lengths:

    repSTRdiff("runparams.txt", [20000:10000:300000...])

Examine the results in runparams.txt.out.final. 

-----------------------------------------------------------

"""
function repSTRdiff(infile::String, hapsizes::Array; keepall=false, showhaps=false, min_uss_score=10)

    sites = []
    n = splitext(infile)
    fn = n[1]

    global s_haps = showhaps
    
    open(infile, "r") do file   #get all loci to run
        for i in eachline(file)

            if occursin(r"^#", i)
                continue
            end

            s = split(i, '\t')
            l = s[4]
            push!(sites, l)

        end
    end


    for i in hapsizes
        getSTRdiff(infile, i, showhaps=s_haps)
    end

    for i in sites
        run(`bash -c "grep $i $fn.*.out |sed 's/^.*\.[0-9]*\.out://' > $fn.$i.all.out" `)
    end

    header = "#ped\tiid\tloc\ttransmitting_gp\ttransmitting_gp_sex\ttransmitting_gp_str\ttransmitting_gp_str_size\tparent_str\tparent_str_size\tstr_allele_change\tsize_change_bp\tphased_offspring_with_de_novo\tfour_gp_match_probabilities\tuniq_soln_score\thapsize\tstr_type"

    OUT = open("$fn.out.final", "w")

    println(OUT, header)

    for i in sites
        best = ""
        open("$fn.$i.all.out") do file

            uss = 0; off = 0;

            for i in eachline(file)

                l = []
                l = split(i, '\t')
                f = String(l[4])
                global ped = l[1]
                global iid = l[2]
                global loc = l[3]

                if f == "na"
                    continue
                end

                if length(l) ≥ 14

                    t = parse(Int, l[14])
                    o = parse(Int, l[12])
                    loc = l[3]

                    if (t + o) > (uss + off)
                        best = i
                        uss = t
                        off = o
                    end

                else

                    best = "$ped\t$iid\t$loc\toptimal solution not found"

                end

            end

            if length(best) < 1
                best = "$ped\t$iid\t$loc\toptimal solution not found"
            end


            function uss_check(best, buss, min_uss_score)
                if (buss ≤ min_uss_score) && (length(best) > 70) && (!occursin(r"\t-9\t", best))
                    b = split(best, '\t')
                    s = parse(Int64, b[14])
                    c = replace(b[13], "[" => "")
                    c = replace(c, "]" => "")
                    p = parse.(Float64, split(c, r", "))
                    d1 = abs(p[1] - p[3])
                    d2 = abs(p[1] - p[4])
                    d3 = abs(p[2] - p[3])
                    d4 = abs(p[2] - p[4])
#                    println("$d1  $d2  $d3  $d4")
                    uss2 = floor(Int64,(minimum([d1, d2, d3, d4]) * 100))
                    b[14] = join([ s, ",", "(", uss2, ")" ],"")
                    best2 = ( join(b, "\t") )
                else
                    best2 = best
                end
                return best2
            end

            if length(best) > 70             #Retest grandparental chromosomes together
                g = split(best, '\t')
                buss = parse(Int64, g[14])
                best = uss_check(best, buss, min_uss_score)
            end

            println(OUT, best)

        end
    end

    close(OUT)

    function clean_dir()
        rt1 = run(`bash -c "rm -f *.out"`);
        rt2 = run(`bash -c "rm -f *_fam* *_trio*"`);
        rt3 = run(`bash -c "rm -f log* error*"`);
        println()
    end

    if keepall == false
        clean_dir()
    elseif keepall == true
        println("INFO: intermediate files for STR sites and haplotype lengths retained.")
    else
        clean_dir()
    end


end
