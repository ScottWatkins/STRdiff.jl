#!/usr/bin/env juila
"""
Run bgzip on a file. The script will either decompress or
compress a file depending on whether infile ends in ".gz".
Returns the name of the outfile.

bgzip_file(infile)

"""
function bgzip_file(infile::String)
    if occursin(r".gz$", infile)
        run(`bash -c "bgzip -d -@8 $infile"`)
        filename = splitext(infile)
        outfile = filename[1]
    elseif ~occursin(r".gz$", infile)
        run(`bash -c "bgzip -@8 $infile"`)
        filename = join([infile, "gz"], ".")
        outfile = filename
    else
        println("ERROR: Could not bgzip or unzip file $filein.\n\n")
        exit()
    end
    return outfile
end
