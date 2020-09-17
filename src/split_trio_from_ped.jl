"""
    (newids, removedids) = split_trio_from_ped("pedigree", "id")

Take input id and remove that id and any parents of that id
from a pedigree file.  Write a cleaned pedigree file and a
pedigree file for the removed ids to disk.  Return strings
containing the ids for the new files.
"""
function split_trio_from_ped(pedigree, id)

    newids = String[]
    dtrioids = String[]
    ped = String[]
    badids = Dict{String,String}()

    open(pedigree) do file

        replace_parent = String[]

        for i in eachline(file)
            l = split(i, '\t')
            if l[2] == id
                badids[l[2]] = "bad"
                badids[l[3]] = "bad"
                badids[l[4]] = "bad"
                replace_parent = [l[1], l[2], "0", "0", l[5], l[6]]
            end
            push!(ped, i)
        end

        OUT = open("$pedigree.fam", "w")
        OUT2 = open("$pedigree.trio", "w")

        for i in ped
            l = split(i, '\t')
            x = l[2]
            if haskey(badids, x)
                push!(dtrioids, x)
                println(OUT2, i)
            else
                push!(newids, x)
                println(OUT, i)
            end
        end

        println(OUT, join(replace_parent, '\t'))
        push!(newids, id)

        close(OUT)
        close(OUT2)

    end

    newids = join(newids, ",")
    dtrioids = join(dtrioids, ",")
    return (newids, dtrioids)

end
