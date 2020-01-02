#ENV["PYTHON"] = "/home/matthewmcallisterdavis/miniconda3/envs/psi/bin/python"
#import Pkg
#Pkg.build("PyCall")
using PyCall
using JuES.Wavefunction
using JuES.CISingles
psi4 = pyimport("psi4")
psi4.core.be_quiet()

mol = psi4.geometry("""
					O
					H 1 1.1
					H 1 1.1 2 104.0
		    symmetry c1
		    """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type" => "pk"))
e,wfn = psi4.energy("hf",mol=mol,return_wfn=true)
Wfn = PyToJl(wfn,Float64,false)
print(do_CIS(Wfn,2,"diag",true))
print(do_CIS(Wfn,2,"svd",true))
print(do_CIS(Wfn,2,"iter",true))
