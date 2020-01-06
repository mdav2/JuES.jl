using Test
using PyCall
#using BenchmarkTools
using JuES.Wavefunction
using JuES.MollerPlesset
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
psi4.set_options(Dict("D_CONVERGENCE" => 10))
mol2 = psi4.geometry("""
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
e2,wfn2 = psi4.energy("hf/sto-3g",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2,Float64,true,false)
mol3 = psi4.geometry("""
					 1 2
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
psi4.set_options(Dict("reference" => "uhf"))
e3,wfn3 = psi4.energy("hf/sto-3g",mol=mol3,return_wfn=true)
JuWfn3 = Wfn(wfn3,Float64,true,false)
@testset "MP2" begin
	@test do_rmp2(JuWfn2) ≈ -0.04915088620016106
	@test do_ump2(JuWfn3) ≈ -0.03588729625230
	@test do_rmp2(JuWfn2) ≈ do_ump2(JuWfn2)
end
