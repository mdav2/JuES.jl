module FCI
using Combinatorics
using LinearAlgebra
using TensorOperations
using JuES
using JuES.Output
using JuES.Wavefunction
using JuES.ConfigurationInteraction
using JuES.ConfigurationInteraction.DetOperations
using JuES.ConfigurationInteraction.MatrixElement

export do_fci

function do_fci(wfn::Wfn; kwargs...)

    # Print intro
    JuES.ConfigurationInteraction.print_header()
    @output "\n    • Computing FCI with the FCI module.\n\n"
    
    # Process options
    for arg in keys(JuES.ConfigurationInteraction.defaults)
        if arg in keys(kwargs)
            @eval $arg = $(kwargs[arg])
        else
            @eval $arg = $(JuES.ConfigurationInteraction.defaults[arg])
        end
    end

    # Integral transformation
    C = wfn.Ca
    gao = wfn.ao_eri
    hao = wfn.hao
    @tensoropt begin
        V[i,a,j,b] := C[σ,b]*C[λ,j]*C[ν,a]*C[μ,i]*gao[μ,ν,λ,σ]
        h[p,q] := C[μ,p]*C[ν,q]*hao[μ,ν] 
    end

    # Check if the number of electrons is even
    nα = wfn.nalpha
    nβ = wfn.nbeta
    nmo = wfn.nmo

    # Create a reference
    zeroth = vcat(repeat([1], nα),repeat([0], nmo-nα))

    perms = permutations(zeroth)

    dets = []

    for αstring in perms
        for βstring in perms
            _det = Determinant(join(αstring), join(βstring))
            push!(dets, _det)
        end
    end

    Ndets = length(dets)

    H = zeros(Ndets, Ndets)
    for (i, D1) in enumerate(dets)
        for (j, D2) in enumerate(dets)
            if excitation_level(D1,D2) > 2
                nothing
            elseif excitation_level(D1,D2) == 2
                H[i,j] = Hd2(D1, D2, V)
            elseif excitation_level(D1,D2) == 1
                αindex = αindex(D1)
                βindex = βindex(D1)
                H[i,j] = Hd1(αindex, βindex, D1, D2, h, V)
            else
                αindex = αindex(D1)
                βindex = βindex(D1)
                H[i,j] = Hd0(αindex, βindex, D1, D2, h, V)
            end
        end
    end

    evals = eigvals(H)

    @output "\n Final FCI Energy: {:15.10f}" evals[1]
                 
end

end #module
