using Test
using PyCall
#using BenchmarkTools
using JuES.Wavefunction
#using JuES.MollerPlesset
using JuES.Transformation
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
JuWfn3 = Wfn(wfn3,Float64,true,false,"uhf")
@testset "Integral Transformation" begin
	@testset "SmokeRHF" begin
		disk = disk_tei_transform(JuWfn2.uvsr,JuWfn2.Ca,"testdisk")
		mem = mem_tei_transform(JuWfn2.uvsr,JuWfn2.Ca)
		@test disk[:,:,:,:] == mem[:,:,:,:]
	end
	@testset "UHF" begin
		disk_tei_transform(JuWfn3.uvsr,JuWfn3.Ca,JuWfn3.Cb,JuWfn3.Ca,JuWfn3.Cb,"TestTEImixed")
		rhf = disk_tei_transform(JuWfn2.uvsr,JuWfn2.Ca,"TestTEIrhf")
		uhf = disk_tei_transform(JuWfn2.uvsr,JuWfn2.Ca,JuWfn2.Ca,JuWfn2.Ca,JuWfn2.Ca,"TestTEIuhf")
		@test rhf[:,:,:,:] ≈ uhf[:,:,:,:]
	end
end
