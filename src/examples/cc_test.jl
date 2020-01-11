using PyCall
using Profile
using InteractiveUtils
using JuES.Wavefunction
using JuES.CoupledCluster
psi4 = pyimport("psi4")
psi4.core.be_quiet()

mol = psi4.geometry("""
		    O
		    H 1 1.1
		    H 1 1.1 2 104.0
		    symmetry c1
		    """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type" => "pk",
					  "d_convergence" => 14))
e,wfn = psi4.energy("hf",mol=mol,return_wfn=true)
refWfn = PyToJl(wfn,Float64,false)
print(do_rccd(refWfn,5))
print(@time do_rccd(refWfn,40))
