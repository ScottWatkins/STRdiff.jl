"""
    get_pedfile_trios(pedfile.txt)

Read a standard format pedigree file and extract
complete trios and return a dictionary of pedigrees
with a list of trios (child,father,mother) in
each pedigree.
"""
function get_pedfile_trios(pedfile)

ptrios = Dict{String, Array}()

open(pedfile) do file

    for i in eachline(file)

        if occursin(r"^[#|P|p]", i)
            continue
        end

        line = split(i, r"\t|\,|\s+")

        if (line[3] != "0") && (line[4] != "0")

            t = join([ line[2], line[3], line[4] ], ",")

            if !haskey(ptrios, line[1])
                ptrios[ line[1] ] = [t]
            else
                push!(ptrios[ line[1] ], t)
            end

        end
    end
end
    peds = collect(keys(ptrios))
    return(ptrios, peds)
end
