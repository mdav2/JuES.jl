module ecRCCSD
using JuES
using JuES.Output
using JuES.Wavefunction

export do_ecrccsd

function do_ecrccsd(wfn::Wfn, detciout::String; kwargs...)

    # Print intro
    JuES.CoupledCluster.print_header()
    @output "\n    • Computing Externally Corrected CCSD with the ecRCCSD module.\n\n"
    
    # Process options
    for arg in keys(JuES.CoupledCluster.defaults)
        if arg in keys(kwargs)
            @eval $arg = $(kwargs[arg])
        else
            @eval $arg = $(JuES.CoupledCluster.defaults[arg])
        end
    end

    # Check if the number of electrons is even
    nelec = wfn.nalpha + wfn.nbeta
    nelec % 2 == 0 ? nothing : error("Number of electrons must be even for RHF. Given $nelec")
    nmo = wfn.nmo
    ndocc = Int(nelec/2)
    nvir = nmo - ndocc

    # Create HF determinant
    ref = [repeat([1], ndocc); repeat([0], nvir)]

    # Get T1 CAS

end

function extract_t1(fname::String; max_exc::Int=4)
    open(fname) do data
        for ln in eachline(data)
            m = match(r"\s*\*\d*\s+([-]?\d+\.\d+)\s+\(.+,.+\)\s+(.+)", ln)
        end
    end
end
    
end #Module
