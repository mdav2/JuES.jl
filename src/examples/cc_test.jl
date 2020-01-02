using PyCall
using Profile
using InteractiveUtils
#include("CoupledCluster.jl")
using CoupledCluster
using Wavefunction
psi4 = pyimport("psi4")
#include("Crutch.jl")
psi4.core.be_quiet()

mol = psi4.geometry("""
		    O
		    H 1 1.1
		    H 1 1.1 2 104.0
		    symmetry c1
		    """)
#mol = psi4.geometry("""
#                    pubchem:ethane
#		    symmetry c1
#		    """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type" => "pk",
					  "d_convergence" => 14))
wfn = init(mol)
println(wfn)
refWfn = PyToJl(wfn,Float64,false)
print(do_rccd(refWfn,5))
print(@time do_rccd(refWfn,40))
#Profile.print()
