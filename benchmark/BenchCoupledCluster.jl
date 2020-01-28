using Test
using PyCall
using BenchmarkTools
using JuES.Wavefunction
using JuES.CoupledCluster: RCCSD,RCCD

psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
psi4.set_num_threads(6)
# > setup
tol = 1E-14
psi4.set_options(Dict("D_CONVERGENCE" => 14,
                      "E_CONVERGENCE" => 12,
					  "scf_type" => "pk"))
mol2 = psi4.geometry("""
                     O
                     H 1 1.1
                     H 1 1.1 2 104.0
					 symmetry c1
					 """)
e2,wfn2 = psi4.energy("hf/cc-pvdz",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2,Float64,false,false)
#println(@btime RCCD.do_rccd(JuWfn2,40,doprint=true))
println(@benchmark RCCSD.do_rccsd(JuWfn2,40,doprint=true))
println(psi4.energy("ccsd/sto-3g",mol=mol2) - e2)
