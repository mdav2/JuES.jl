using Test
using PyCall
using JuES.Wavefunction
using JuES.CoupledCluster
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
tol = 1E-14
mol2 = psi4.geometry("""
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type"=>"pk","d_convergence"=>14))
e,wfn2 = psi4.energy("hf/sto-3g",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2)
@testset "CoupledCluster" begin
	@testset "Smoke" begin
		@test do_rccd(JuWfn2,40) - -0.07015050066089029 < tol
	end
end
