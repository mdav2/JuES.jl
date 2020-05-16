#run(`julia mk.jl`)
#run(`julia aux.jl`)
#addprocs(1)
#using JLD
printdo = false
using Distributed
addprocs(1)
using JuES
@everywhere using JuES.Wavefunction
using JuES.CoupledCluster: RCCSD,RCCD,DFRCCD,mRCCD,mRCCSD

using LinearAlgebra
BLAS.set_num_threads(6)
# > setup
tol = 1E-14
@everywhere function psi_setup()
    include("Psi4.jl")
    Psi4.psi4.core.be_quiet() #turn off output
    Psi4.psi4.set_num_threads(6)
    Psi4.psi4.set_options(Dict("D_CONVERGENCE" => 14,
                          "E_CONVERGENCE" => 12,
    					  "scf_type" => "pk"))
    mol2 = Psi4.psi4.geometry("""
                         O
                         H 1 1.1
                         H 1 1.1 2 104.0
    					 symmetry c1
    					 """)
    e2,wfn2 = Psi4.psi4.energy("hf/cc-pvtz",mol=mol2,return_wfn=true)
    mints = Psi4.psi4.core.MintsHelper(wfn2.basisset())
    JuWfn2 = Wfn(wfn2,mints)
    return JuWfn2
end
JuWfn2 = psi_setup()
println(@btime RCCSD.do_rccsd(JuWfn2; doprint=printdo))
println(@btime mRCCSD.do_rccsd(JuWfn2; doprint=printdo))
