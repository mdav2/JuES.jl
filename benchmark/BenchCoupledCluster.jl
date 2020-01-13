using Test
using PyCall
using BenchmarkTools
using JuES.Wavefunction
using JuES.CoupledCluster

psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
psi4.set_options(Dict("D_CONVERGENCE" => 10,
					  "scf_type" => "pk"))
mol2 = psi4.geometry("""
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
e2,wfn2 = psi4.energy("hf/cc-pvdz",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2,Float64,false,false)
println(@benchmark do_rccd(JuWfn2))
