using Test
using PyCall
#using BenchmarkTools
using JuES.Wavefunction
using JuES.MollerPlesset
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
mol2 = psi4.geometry("""
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
e2,wfn2 = psi4.energy("hf/sto-3g",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2)
mol3 = psi4.geometry("""
					 1 2
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
psi4.set_options(Dict("reference" => "uhf"))
e3,wfn3 = psi4.energy("hf/sto-3g",mol=mol3,return_wfn=true)
JuWfn3 = Wfn(wfn3)
@testset "Integral Transformation" begin
	@testset "SmokeRHF" begin
		disk = transform_tei(JuWfn2.uvsr,JuWfn2.Ca)
		mem = transform_tei2(JuWfn2.uvsr,JuWfn2.Ca)
		@test disk[:,:,:,:] == mem[:,:,:,:]
	end
	@testset "UHF" begin
		transform_tei(JuWfn3.uvsr,JuWfn3.Ca,JuWfn3.Cb,JuWfn3.Ca,JuWfn3.Cb)
		rhf = transform_tei(JuWfn2.uvsr,JuWfn2.Ca)
		uhf = transform_tei(JuWfn2.uvsr,JuWfn2.Ca,JuWfn2.Ca,JuWfn2.Ca,JuWfn2.Ca)
		@test rhf[:,:,:,:] == uhf[:,:,:,:]
	end
end
