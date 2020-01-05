using Test
using PyCall
using JuES.Wavefunction
using JuES.CISingles
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
tol = 1E-14
mol = psi4.geometry("""
					H
					H 1 1.0
					symmetry c1
					""")
psi4.set_options(Dict("scf_type" => "pk","d_convergence"=>14))
e,wfn = psi4.energy("hf/sto-3g",mol=mol,return_wfn=true)
JuWfn = PyToJl(wfn,Float64,false)
@testset "CISingles" begin
	@testset "Smoke" begin
		@test do_CIS(JuWfn,1)[1] - 0.320236855893771 < tol 
	end
end
