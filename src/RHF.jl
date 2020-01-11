struct RHFWfn
    molecule::PyObject
    basis::PyObject
    mints::PyObject
    C::Array{Float64,2}
    h::Array{Float64,2}
    S::Array{Float64,2}
    D::Array{Float64,2}
    vnuc::Float64
end

function RHFWfn(molecule::PyObject,basis::String="STO-3G")
    vnuc = molecule.nuclear_repulsion_energy()
    basis = psi4.core.BasisSet.build(molecule,fitrole="ORBITAL",target="STO-3G")
    mints = psi4.core.MintsHelper(basis)
    S = mints.ao_overlap().np
    h = mints.ao_kinetic().np + mints.ao_potential().np
    C = zeros(basis.nao(),basis.nao())
    D = zeros(basis.nao(),basis.nao())
    RHFWfn(molecule,basis,mints,C,h,S,D,vnuc)
end

function RHFCompute(wfn::RHFWfn)
    ao_eri = RHFWfn.ao_eri().np

end
"""
    RHFEnergy
computes the RHF energy given a (not necessarily converged) core and
two electron integrals in the AO basis, as well as the desired density matrix.
## paramters
    h::Array{Float64,2} -> core hamiltonian array in AO basis
    g::Array{Float64,4} -> two electron integrals (eq 4b)
                        ->      g[μ,ν] = Σ(ρ,σ) ⟨μρ|νσ⟩ - ⟨μρ|σν⟩
    D::Array{Float64,2} -> density matrix (eq 2b) second C is actually Cconj
                        ->      D[μ,ν] = 2Σ(i;1:nocc) C[μ,i]*C[ν,i]
## outputs
    E::Float            -> SCF energy
"""
function RHFEnergy(h,g,D)
    nao = size(h)[1]
    R = 1:nao
    v = zeros(size(D))
    #form v (eqn 4b)
    for μ in R
        for ν in R
            temp = 0.0
            for ρ in R
                for σ in R
                    temp += (g[μ,ν,ρ,σ] - (1/2)*g[μ,σ,ρ,ν])*D[σ,ρ]
                end
            end
            v[μ,ν] = temp
        end
    end
    #compute the energy
    E = 0.0
    for ν in R
        for μ in R
            E += (h[μ,ν] + (1/2)*v[μ,ν])*D[ν,μ]
        end
    end
    return E
end

"""
    form_D

form the density matrix from an AO->MO coefficient matrix.
"""
function form_D(nocc,C)
    nao = size(C)[1]
    R = 1:nao
    D = zeros(nao,nao)
    for μ in R
        for ν in R
        end
    end
end
