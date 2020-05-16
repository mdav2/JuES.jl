using JuES
using JLD

include("Psi4.jl")

function psi_setup()
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
    JuWfn2 = JuES.Wavefunction.Wfn{Float64}(wfn2,mints)
    return JuWfn2
end

wfn = psi_setup()
save("/tmp/mywfn.jld","mywfn",wfn)
