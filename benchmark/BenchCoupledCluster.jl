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
					  "scf_type" => "pk"))
mol2 = psi4.geometry("""
                     pubchem:butane
					 symmetry c1
					 """)
e2,wfn2 = psi4.energy("hf/cc-pvdz",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2,Float64,false,false)
println(@btime RCCD.do_rccd(JuWfn2,40,doprint=true))
println(@btime RCCSD.do_rccsd(JuWfn2,40,doprint=true))
