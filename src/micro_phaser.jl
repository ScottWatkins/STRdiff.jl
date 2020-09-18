"""
    micro_phaser(mypedigree.vcf.gz, mypedfile.ped)

Run Beagle 4.0 on a vcf file using a pedigree to inform
the phasing.  A map file is not used, so the region should
be small, in the 1-150 kb range.

"""
function micro_phaser(vcf, ped)

    if Sys.which("beagle.r1399.jar") == nothing
        error("beagle.r1399.jar is required in your path")
    else
        beagle = Sys.which("beagle.r1399.jar")
    end

    printstyled("TEST1", color=:red)

    run(`bash -c "java -Xmx1g -jar $beagle gt=$vcf ped=$ped out=phased 1>log 2>error" `)
    run(`mv phased.vcf.gz $vcf`)
    run(`rm phased.log phased.warnings`)  #optional directory cleaning

end
