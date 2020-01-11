using PyCall
using Wavefunction
using Direct
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
myWfn = DirectWfn(wfn)
println(ao(myWfn.mints,myWfn.basis,1,1,1,1))
out = func_shell(myWfn.basis)
println(out)
