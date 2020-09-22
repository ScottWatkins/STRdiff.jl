"""
    repSTRdiff("infile.txt", [hapsize_1, hapsize_2, ... hapsize_N];
        keepall=false, showhaps=false, min_uss_score=10)

Run getSTRdiff for multiple user specified haplotype sizes.
Retrieve the best estimates from all runs and output a final table.
See getSTRdiff for additional details.

Options:

    keepall        Keep all intermediate files [false]

    showhaps       Show the offspring and de novo trio haplotypes [false]

    min_uss_score  Minimum value for the unique solution score.  For values
                   below this number, both chromosomes from each of the
                   grandparents will be tested and the minimum uss score for
                   the comparisons will be shown parenthetically, thus allowing
                   a potential transmitting parent to be found even though
                   no single chromosome solution produced a high score [10].

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

            function uss_check(best)

                if length(best) > 70

                    b = split(best, '\t')
                    s = parse(Int64, b[14])

                    if s ≤ min_uss_score
                        c = replace(b[13], "[" => "")
                        c = replace(c, "]" => "")
                        p = parse.(Float64, split(c, r", "))
                        d1 = abs(p[1] - p[3])
                        d2 = abs(p[1] - p[4])
                        d3 = abs(p[2] - p[3])
                        d4 = abs(p[2] - p[4])
                        if (d1 ≥ 0.1 && d2 ≥ 0.10 && d3 ≥ 0.1 && d4 ≥ 0.1)
                            uss2 = floor(Int64,(minimum([d1, d2, d3, d4]) * 100))
                        end
                        
                        b[14] = join([ s, ",", "(", uss2, ")" ],"")
                        best = ( join(b, "\t") )

                    end

                end

                return best

            end

            best = uss_check(best)
            
            println(OUT, best)

        end
    end

    close(OUT)

    function clean_dir()
        rt1 = run(`bash -c "rm -f *.out"`);
        rt2 = run(`bash -c "rm -f *_fam* *_trio*"`);
        rt3 = run(`bash -c "rm -f log*"`);
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
