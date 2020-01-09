using Test
using PyCall
using BenchmarkTools
using JuES.Wavefunction
using JuES.MollerPlesset

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
JuWfn2 = Wfn(wfn2,Float64,true,false)
mol3 = psi4.geometry("""
					 1 2
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
psi4.set_options(Dict("reference" => "uhf"))
e3,wfn3 = psi4.energy("hf/cc-pvdz",mol=mol3,return_wfn=true)
println(@benchmark Wfn(wfn3,Float64,true,true))
println(@benchmark Wfn(wfn3,Float64,true,false))
JuWfn3 = Wfn(wfn3,Float64,true,true)
JuWfn4 = Wfn(wfn3,Float64,true,false)
println(@benchmark do_rmp2(JuWfn2))
println(@benchmark do_ump2(JuWfn3))
println(@benchmark do_ump2(JuWfn4))
