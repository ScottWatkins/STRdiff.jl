"""
    check_vcf_sites("infile.vcf.gz", "datafile.txt"; c=1)

Check that variant sites (e.g. 9:43985033) listed in a file
are found in a given vcf file. List the sites not found in the
vcf file.  Write a new file that contains only the sites found
in the vcf file.

Options:

c=n    Specify the column number n that contains the sites [1].

Dependencies: DelimitedFiles
"""
function check_vcf_sites(infile::String, datafile::String; c=1)

    print("Reading all sites from $datafile ... ")
    filein = readdlm(datafile, '\t')
    filesites = Set( readdlm(datafile, '\t')[:,c] )
    println("(e.g. ", rand(filesites), ")" , " ... done.")

    print("Reading all sites from $infile ... ")
    vcfsites = Set( split( read(`bcftools query -f '%CHROM:%POS,' $infile`, String), ',' ) )
    println("done.")

    fileonly = setdiff(filesites, vcfsites)
    isect = intersect(filesites, vcfsites)

    ni = length(isect)
    n = length(fileonly)

    println("Found $ni sites in vcf file.")
    println("Found $n sites not in vcf file:")
    println(fileonly)

    gsites = filein[ findall(in(isect), filein[:,c]), : ]

    writedlm("$datafile.filtered", gsites, '\t')
    println("Wrote filtered file to $datafile.filtered.")

end
