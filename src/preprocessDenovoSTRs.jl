"""
    preprocessDenovoSTRs("infile.tsv"; window=100000)

A script to extract a region around a specified STR site
where a event of interest (e.g. mutation) may have occurred.

Options:

path    The directory contaning the master vcf file(s) and
        and pedigree file(s) [default: current dir].
        The pedigree files must be named myped.ped (e.g. 1463.ped)

window  The number of bp ± each side of the STR site to construct the
        haplotype [default: 100000 bp].

Returns a list containing the extracted vcf, pedfile name,
and STR location for input into the STR_iht.jl script and a small
phased vcf of the extracted regions.
"""
function preprocessDenovoSTRs(infile; kwargs...)

outdat = []

if Sys.which("bcftools")  == nothing
    error("bcftools is required. Please add bcftools to your path.")
end

kw = Dict(kwargs)

#get path (or set current path to) directory with phased vcfs
path = get(kw, :path, pwd()*"/")

#set haplotype range (each side)
window = get(kw, :window, 100000)


        i = infile
        a = split(i, '\t')
        ped = a[1]
        vcf = a[2]
        pname = splitext(a[1])
        iid = a[3]
        site = a[4]
        (chr, pos) = split(a[4], ':')
        pos = parse(Int64, pos)
        start = pos - window
        stop  = pos + window
        range = join([chr, start, stop], ":", "-")

        (famids, tids) = split_trio_from_ped(ped, iid)

        f_outfile = join([pname[1], "_", "fam", "_", a[4], ".vcf.gz"], "")
        t_outfile = join([pname[1], "_", "trio", "_", a[4], ".vcf.gz"], "")

        if !isfile("$vcf.csi")
            println("INFO: indexing the vcf file $vcf")
            run(`bcftools index $vcf`)
        end

        rs_f = success(run(`bcftools view -s $famids -r $range --min-ac 3 -e'ALT="*" || QUAL < 100000 || FORMAT/DP < 10 || MQ < 30 || GQ < 30' $vcf -Oz -o $f_outfile`, wait=true))

        if rs_f == true
            run(`bcftools index -f $f_outfile`, wait=true)
        end

        run(pipeline(`bcftools view -H $f_outfile`, `cut -f1,2`, "famsites.txt"))

        t_tmp = success(run(`bcftools view -s $tids -R famsites.txt $vcf -Oz -o "triovcftmp.vcf.gz"`, wait=true) )

        t_tmp2 = success(run(`bcftools norm -d none triovcftmp.vcf.gz -Oz -o "triovcftmp2.vcf.gz"`))

        rs_t = success(run(`bcftools sort triovcftmp2.vcf.gz -Oz -o $t_outfile`))

        if rs_t == true
            run(`bcftools index -f $t_outfile`, wait=true)
        end

        rs_s = success(run(`bcftools view -s $tids -r $site $t_outfile -Oz -o $site.only.vcf.gz`, wait=true) )

        if rs_s == true
            run(`bcftools index -f $site.only.vcf.gz`, wait=true)
        end

        run(`bcftools filter -i 'FMT/DP > 99999999' -S . -r $site $site.only.vcf.gz -Oz -o missing.vcf.gz`, wait=true)

        function getSTRmut(iid)
            child = split(read(`bcftools view -H -s $iid $site.only.vcf.gz `, String), '\t')
            cs = Set([child[end][1], child[end][3]])
            cs_geno = join([ child[end][1], child[end][3] ], "/")
            par = split(read(`bcftools view -H -s ^$iid $site.only.vcf.gz `, String), '\t')
            ps_ids = split(chomp(read(`bash -c "bcftools view -h -s ^$iid $site.only.vcf.gz | tail -n1 | cut -f10-"`, String)), '\t')
            ps = Set( [ par[end][1], par[end][3], par[end-1][1], par[end-1][3] ])

            p1 = ps_ids[1]
            p2 = ps_ids[2]
            par_genos = (p1, par[end-1][1], par[end-1][3], p2, par[end][1], par[end][3] )

            mutSTR = collect(setdiff(cs, ps))

            raw_genos = join([ p1, ": ", par[end-1][1], "/", par[end-1][3], ": ", p2, ": ", par[end][1], "/", par[end][3], " child: ", cs_geno], "")

            if length(mutSTR) == 0
                error("no unique STR allele in child: $par_genos -> $cs_geno")
            end

            mutSTR = parse(Int, mutSTR[1])
            return mutSTR, par_genos, raw_genos
        end

        (mutSTR, par_genos, raw_genos) = getSTRmut(iid)   #child STR allele ∉ parents

        rs_ms = success(run(`cp missing.vcf.gz $site.only.vcf.gz`, wait=true) )

        if rs_ms == true
            run(`rm -f missing.vcf.gz`)
            run(`bcftools index -f $site.only.vcf.gz`)
        end

        function micro_phaser(vcf, ped)

            if Sys.which("beagle.r1399.jar") == nothing
                error("beagle.r1399.jar is required in your path")
            else
                beagle = Sys.which("beagle.r1399.jar")
            end

            run(`java -Xmx1g -jar $beagle gt=$vcf ped=$ped out=phased`)
            run(`bcftools index -f phased.vcf.gz`)
            run(`mv phased.vcf.gz $vcf`)
            run(`mv phased.vcf.gz.csi $vcf.csi`)

        end

        famped = join([ped, "fam"], '.')
        trioped = join([ped, "trio"], '.')

        micro_phaser(f_outfile, famped)
        micro_phaser(t_outfile, trioped)

        tm = success(run(`bcftools concat $site.only.vcf.gz $t_outfile -Oz -o trio_merged.vcf.gz`) )

        tms = success(run(`bcftools sort trio_merged.vcf.gz -Oz -o trio_merged_sort.vcf.gz`) )

        if tms == true
            run(`bcftools index  -f trio_merged_sort.vcf.gz`)
        end

        run(`bash -c "mv trio_merged_sort.vcf.gz $t_outfile"`)
        run(`bcftools index -f $t_outfile`)

        #Need to rematch  and renorm the sites due to misinheritances
        #thrown out when phasing in some trios
        run(pipeline(`bcftools view -H $t_outfile`, `cut -f1,2`, "famsites_t.txt"))

        f_tmp = success(run(`bcftools view -s $famids -R famsites_t.txt $f_outfile -Oz -o fam_rematchedtmp.vcf.gz`, wait=true) )

        f_tmp2 = success(run(`bcftools norm -d none fam_rematchedtmp.vcf.gz -Oz -o fam_rematchedtmp2.vcf.gz`))

        rs_ft = success(run(`bcftools sort fam_rematchedtmp2.vcf.gz -Oz -o $f_outfile`))

        if rs_ft == true
            run(`bcftools index -f $f_outfile`)
        end

        run(`bash -c "rm -f trio_merged.vcf.gz $site.only.vcf.gz trio_merged_sort.vcf.gz.csi  $site.only.vcf.gz.csi famsites.txt phased.* famsites_t.txt *tmp*" `)

        printstyled("INFO: unphased de novo trio genotypes: ", raw_genos, color=:green)

        STR_location = a[4]

        push!(outdat, [f_outfile, famped, STR_location, iid, t_outfile, trioped, mutSTR, par_genos] )

        flines = strip(read(`bash -c "bcftools view -H $f_outfile | wc -l" `, String))
        tlines = strip(read(`bash -c "bcftools view -H $f_outfile | wc -l" `, String))

        if flines != tlines
            error("Varinats in vcfs not equal for $STR_location, got fam vcf=$flines, trio vcf=$tlines. Suggest removing $STR_location.")
            exit()
        end

    return outdat

end
