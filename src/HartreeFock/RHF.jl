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
    grad::Bool
    hess::Bool
    GradN::Array{Float64,2}
    GradS::Array{Float64,2}
    GradSp::Array{Float64,2}
    GradV::Array{Float64,2}
    GradT::Array{Float64,2}
    GradJ::Array{Float64,2}
    GradK::Array{Float64,2}
    Grad::Array{Float64,2}
end

function RHFWfn(molecule::PyObject,basis::String="STO-3G";debug=false,grad=false,hess=false)
    dummy2 = Array{Float64}(undef,0,0)
    dummy4 = Array{Float64}(undef,0,0,0,0)
    natoms = molecule.natom()
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
    if ! grad
        GradN = dummy2
        GradS = dummy2
        GradSp = dummy2
        GradV = dummy2
        GradT = dummy2
        GradJ = dummy2
        GradK = dummy2
        Grad = dummy2
        
    end
    RHFWfn(molecule,basis,mints,C,H,S,A,D,vnuc,nelec/2,grad,hess,GradN,GradS,
          GradSp,GradV,GradT,GradJ,GradK,Grad)
end

function RHFCompute(wfn::RHFWfn;doprint=false,maxit=50,Etol=1E-7,Dtol=1E-7)
    print_header()
    @output "    executing RHF\n"
    @output "    computing AO basis integrals ... "
    t = @elapsed I = wfn.mints.ao_eri().np
    G = 2*I - permutedims(I,[1,3,2,4])
    @output "done in {:>5.2f}s\n" t
    @output "    Forming initial Fock matrix ... "
    t = @elapsed begin
        Ft = transpose(wfn.A)*wfn.H*wfn.A
        e,Ct = eigen(Ft)
        C = wfn.A*Ct
        Co = C[:,1:wfn.ndocc]
    end
    @output "done in {:>5.2f}s\n" t
    @tensor D[u,v] := Co[u,m]*Co[v,m]
    @tensor F[m,n] := D[r,s]*G[m,n,r,s]
    F += wfn.H
    E = 0#RHFEnergy(D,wfn.H,F) + wfn.vnuc
    #if doprint println("@RHF 0 $E") end

    @output "\n"
    @output " Iter.   {:<20} {:>11} {:>11} {:>8}\n" "E[RHF]" "dE" "√|D|²" "t"
    @output repeat("~",80)*"\n"
    t = @elapsed for i in 1:maxit
        t_iter = @elapsed begin
            @tensor F[m,n] := wfn.H[m,n] + D[r,s]*G[m,n,r,s]
            Eelec = RHFEnergy(D,wfn.H,F)
            Enew = Eelec + wfn.vnuc
            Ft = transpose(wfn.A)*F*wfn.A
            Ft = Symmetric(Ft)
            e,Ct = eigen(Ft)#,sortby = x->-abs(x))
            C = wfn.A*Ct
            Co = C[:,1:wfn.ndocc]
            @tensor Dnew[u,v] := Co[u,m]*Co[v,m]
            dD = Dnew - D
            Drms = sqrt(sum(dD)^2)
            dE = Enew - E
            D = Dnew
            E = Enew
        end
        #if doprint println("@RHF $i $E $dE $Drms") end
        @output "    {:<3} {:>20.17f} {:>11.3e} {:>11.3e} {:>8.2f}\n" i E dE Drms t_iter
        if (dE < Etol) & (Drms < Dtol)
            break
        end
    end
    @output repeat("~",80)*"\n"
    @output "    RHF done in {:>5.2f}s\n" t
    @output "    @E[RHF] = {:>20.17f}" E
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

