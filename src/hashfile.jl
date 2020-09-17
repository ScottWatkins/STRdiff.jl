"""
    hashfile("inputfile", keycol, valcols)

Create a dictionary from a tsv file with some
column of the file as the lookup key and values
from columns specified in an array. 
kwargs: delim="\t"
"""
function hashfile(file::String, keycolumn::Int, valcols::Array; delim="\t")

    fhash = Dict{Any,Any}()
    d = Regex(delim)
    
    open(file) do f

        for i in eachline(f)
            line = split(i, d)
            fhash[line[keycolumn]] = line[valcols]
        end

    end
    
    return fhash

end
