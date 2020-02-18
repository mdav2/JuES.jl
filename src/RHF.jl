struct RHFWfn
    molecule::PyObject
    basis::PyObject
    mints::PyObject
    C::Array{Float64,2}
    H::Array{Float64,2}
    S::Array{Float64,2}
    A::Array{Float64,2}
    D::Array{Float64,2}
    vnuc::Float64
    ndocc::Int
end

function RHFWfn(molecule::PyObject,basis::String="STO-3G";debug=false)
    vnuc = molecule.nuclear_repulsion_energy()
    basis = psi4.core.BasisSet.build(molecule,fitrole="ORBITAL",target="STO-3G")
    mints = psi4.core.MintsHelper(basis)
    S = mints.ao_overlap().np
    A = mints.ao_overlap()
    A.power(-0.5,1e-16)
    A = A.np
    if debug println("Setup done") end
    T = mints.ao_kinetic().np 
    if debug println("made T") end
    V = mints.ao_potential().np
    if debug println("made V") end
    H = T+V
    if debug println("Made H") end
    C = zeros(basis.nao(),basis.nao())
    if debug println("Made C") end
    D = zeros(basis.nao(),basis.nao())
    if debug println("Made D") end
    nelec = -1*molecule.molecular_charge()
    for i in 0:molecule.natom()-1
        nelec += molecule.charge(i)
    end
    if debug println("Computed nelec") end
    if nelec%2 != 0
        return false
    end
    RHFWfn(molecule,basis,mints,C,H,S,A,D,vnuc,nelec/2)
end

function RHFCompute(wfn::RHFWfn;doprint=false,maxit=50,Etol=1E-7,Dtol=1E-7)
    I = wfn.mints.ao_eri().np
    G = 2*I - permutedims(I,[1,3,2,4])
    Ft = transpose(wfn.A)*wfn.H*wfn.A
    e,Ct = eigen(Ft)
    C = wfn.A*Ct
    Co = C[:,1:wfn.ndocc]
    @tensor begin
        D[u,v] := Co[u,m]*Co[v,m]
    end
    @tensor begin
        F[m,n] := D[r,s]*G[m,n,r,s]
    end
    F += wfn.H
    E = RHFEnergy(D,wfn.H,F) + wfn.vnuc
    if doprint println("@RHF 0 $E") end
    for i in 1:maxit
        @tensor begin
            F[m,n] := wfn.H[m,n] + D[r,s]*G[m,n,r,s]
        end
        Eelec = RHFEnergy(D,wfn.H,F)
        Enew = Eelec + wfn.vnuc
        Ft = transpose(wfn.A)*F*wfn.A
        Ft = Symmetric(Ft)
        e,Ct = eigen(Ft)#,sortby = x->-abs(x))
        C = wfn.A*Ct
        Co = C[:,1:wfn.ndocc]
        @tensor begin
            Dnew[u,v] := Co[u,m]*Co[v,m]
        end
        dD = Dnew - D
        Drms = sqrt(sum(dD)^2)
        dE = Enew - E
        D = Dnew
        E = Enew
        if doprint println("@RHF $i $E $dE $Drms") end
        if (dE < Etol) & (Drms < Dtol)
            break
        end
    end
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
function RHFEnergy(D,H,F)
    temp = H + F
    @tensor begin
        E[] := D[m,n]*temp[m,n]
    end
    return E[]
end

