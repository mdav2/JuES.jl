using PyCall
using BenchmarkTools
using Wavefunction
using MollerPlesset
psi4 = pyimport("psi4")
np = pyimport("numpy")
mol = psi4.geometry("""
					O
				    H 1 1.1
			        H 1 1.1 2 104.0
					symmetry c1""")
psi4.core.be_quiet()

bbasis = "sto-3g"
psi4.set_options(Dict("scf_type" => "pk","basis" => bbasis))
ehf,wfn = psi4.energy("hf",mol=mol,return_wfn=true)
Wfn = PyToJl(wfn,Float64,false)
println(do_rmp2(Wfn))
