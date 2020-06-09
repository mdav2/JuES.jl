module FCI
using Combinatorics
using LinearAlgebra
using SparseArrays
using Arpack
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
            println(arg)
        end
    end

    @output "SCF Energy: {:10.10f}\n" wfn.energy

    nmo = wfn.nmo

    if wfn.nalpha ≠ wfn.nbeta
        error("Open-shell case not implemented yet.")
    end

    act_elec = wfn.nalpha + wfn.nbeta - 2*frozen

    if act_elec ≤ 0
        error("Invalid number of frozen orbitals ($frozen) for $(nα+nβ) electrons.")
    end

    # Active = -1 means FCI, with frozen
    if active == -1
        @eval active = $nmo - $frozen
    end

    if active ≤ act_elec/2
        error("Number of activ orbitals ($active) too small for $(act_elec) activ electrons")
    end

    @output "\n →  ACTIVE SPACE\n"
    @output "Frozen Orbitals:  {:3d}\n" frozen
    @output "Active Electrons: {:3d}\n" act_elec
    @output "Active Orbitals:  {:3d}\n" active
    

    @output "Transforming integrals...\n"

    # Integral transformation
    C = wfn.Ca
    gao = wfn.ao_eri
    hao = wfn.hao
    @tensoropt begin
        V[i,a,j,b] := C[σ,b]*C[λ,j]*C[ν,a]*C[μ,i]*gao[μ,ν,λ,σ]
        h[p,q] := C[μ,p]*C[ν,q]*hao[μ,ν] 
    end

    @output "Done.\n"


    dets = get_determinants(act_elec, active, nmo, frozen)
    Ndets = length(dets)
    @output "Number of Determinants: {:10d}\n" Ndets

    @output "Computing using Full Hamiltonian\n"

    @time begin
        H = zeros(Ndets, Ndets)
        H = get_dense_hamiltonian_matrix(dets, H, h, V)
    end
    @time evals, v = eigs(H, nev=1)
    println(Base.summarysize(H))

    @output "\n Final FCI Energy: {:15.10f}\n" evals[1]+wfn.vnuc
    #@output "Time: {:10.5f}\n" t


    @output "Computing using Sparse Hamiltonian\n"

    @time begin
        Hs = get_sparse_hamiltonian_matrix(dets, h, V)
    end
    @time evals, v = eigs(Hs, nev=1)
    println(Base.summarysize(Hs))

    @output "\n Final FCI Energy: {:15.10f}\n" real(evals[1])+wfn.vnuc
    #@output "Time: {:10.5f}\n" t

end

function get_determinants(Ne::Int, No::Int, nmo::Int, nfrozen::Int)

    Nae = Int(Ne/2)
    occ_string = repeat('1', nfrozen)
    vir_string = repeat('0', nmo-nfrozen-Nae)

    zeroth = repeat('1', Nae)*repeat('0', No-Nae)

    perms = multiset_permutations(zeroth, length(zeroth))

    dets = Determinant[]
    for αstring in perms
        for βstring in perms
            α = occ_string*join(αstring)*vir_string
            β = occ_string*join(βstring)*vir_string
            _det = Determinant(α, β)
            push!(dets, _det)
        end
    end

    # Determinant list is sorted by its excitation level w.r.t the first determinant (normally HF)
    sort!(dets, by=d->excitation_level(dets[1], d))

    return dets
end


function get_dense_hamiltonian_matrix(dets::Array{Determinant,1}, H::AbstractArray{Float64,2}, h::Array{Float64,2}, V::Array{Float64,4})

    Ndets = length(dets)

    for i in 1:length(dets)
        D1 = dets[i]
        αind = αindex(D1)
        βind = βindex(D1)
        for j in i:length(dets)
            D2 = dets[j]
            el = excitation_level(D1,D2)
            if el > 2
                nothing
            elseif el == 2
                H[i,j] = Hd2(D1, D2, V)
            elseif el == 1
                H[i,j] = Hd1(αind, βind, D1, D2, h, V)
            else
                H[i,j] = Hd0(αind, βind, h, V)
            end
        end
    end

    return Symmetric(H)
end

function get_sparse_hamiltonian_matrix(dets::Array{Determinant,1}, h::Array{Float64,2}, V::Array{Float64,4})

    tol = 1e-9

    Ndets = length(dets)

    vals = Float64[]
    ivals = Int64[]
    jvals = Int64[]

    for i in 1:length(dets)
        D1 = dets[i]
        αind = αindex(D1)
        βind = βindex(D1)
        for j in i:length(dets)
            D2 = dets[j]
            el = excitation_level(D1,D2)
            if el > 2
                nothing
            elseif el == 2
                elem = Hd2(D1, D2, V)
                if abs(elem) > tol
                    push!(vals, elem)
                    push!(ivals, i)
                    push!(jvals, j)
                end
            elseif el == 1
                elem = Hd1(αind, βind, D1, D2, h, V)
                if abs(elem) > tol
                    push!(vals, elem)
                    push!(ivals, i)
                    push!(jvals, j)
                end
            else
                elem = Hd0(αind, βind, h, V)
                if abs(elem) > tol
                    push!(vals, elem)
                    push!(ivals, i)
                    push!(jvals, j)
                end
            end
        end
    end

    return Symmetric(sparse(ivals, jvals, vals))
end

function print_matrix(M)

    m,n = size(M)

    for i in 1:m
        println(M[i,:])
    end
end

function sparsity(M::AbstractArray)

    return count(i->abs(i)>1e-10, M)/length(M)

end

end #module
