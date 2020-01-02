#ENV["PYTHON"] = "/home/matthewmcallisterdavis/miniconda3/envs/psi/bin/python"
#import Pkg
#Pkg.build("PyCall")
using PyCall
using Wavefunction
using CISingles
psi4 = pyimport("psi4")
psi4.core.be_quiet()
#include("Crutch.jl")

mol = psi4.geometry("""
					pubchem:propanol
		    symmetry c1
		    """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type" => "pk"))
wfn = init(mol)
Wfn = PyToJl(wfn,Float64,false)
print(do_CIS(Wfn,2,"davidson",true))
print(do_CIS(Wfn,2,"diag",true))
print(do_CIS(Wfn,2,"svd",true))
