"""
    repSTRdiff("infile.txt", [hapsize_1, hapsize_2, ... hapsize_N];
        keepall=false, showhaps=false)

Run getSTRdiff for multiple user specified haplotype sizes.
Retrieve the best estimates from all runs and output a final table.
See getSTRdiff for additional details.

Options:

    keepall    Keep all intermediate files [true|false]
    showhaps   Show the offspring and de novo trio haplotypes [true|false]

"""
function repSTRdiff(infile::String, hapsizes::Array; keepall=false, showhaps=false)

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
