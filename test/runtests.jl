using Test
using PyCall
using JuES.Wavefunction
using JuES.CISingles
using JuES.CoupledCluster
using JuES.DiskTensors
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output
# > setup
tol = 1E-14
mol = psi4.geometry("""
					H
					H 1 1.0
					symmetry c1
					""")
psi4.set_options(Dict("scf_type" => "pk","d_convergence"=>14))
e,wfn = psi4.energy("hf/sto-3g",mol=mol,return_wfn=true)
JuWfn = PyToJl(wfn,Float64,false)

mol2 = psi4.geometry("""
					 O
					 H 1 1.1
					 H 1 1.1 2 104.0
					 symmetry c1
					 """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type"=>"pk","d_convergence"=>14))
e,wfn2 = psi4.energy("hf/sto-3g",mol=mol2,return_wfn=true)
JuWfn2 = PyToJl(wfn2,Float64,false)
include("TestWavefunction.jl")
# < setup

#@testset "Wavefunction" begin
#	@testset "Smoke" begin
#		@testset "Attributes" begin
#			@test JuWfn.nalpha == 1 && JuWfn.nbeta == 1
#		end
#		@testset "Integrals" begin
#			@test (JuWfn.uvsr[1,1,1,1] - 0.7746059439198979) < tol
#			@test (JuWfn.uvsr[2,1,2,2] - 0.3093089669634818) < tol
#			@test (JuWfn.uvsr[1,1,2,2] - 0.4780413730018048) < tol
#		end
#	end
#end
@testset "CISingles" begin
	@testset "Smoke" begin
		@test do_CIS(JuWfn,1)[1] - 0.320236855893771 < tol 
	end
end
@testset "CoupledCluster" begin
	@testset "Smoke" begin
		@test do_rccd(JuWfn2,40) - -0.07015050066089029 < tol
	end
end
include("TestDiskTensors.jl")
#@testset "DiskTensors" begin
#	@testset "Smoke" begin #checking constructors
#		@testset "DiskVector" begin
#			N = 1000
#			x = DiskVector("/tmp/testdv1.h5",Float64,N,"w")
#			y = DiskVector("/tmp/testdv2.h5",Float64,N,"w")
#			@test true
#		end
#		@testset "DiskMatrix" begin
#			N = 1000
#			a = DiskMatrix("/tmp/matrix1.h5",Float64,N,N,"w")
#			b = DiskMatrix("/tmp/matrix2.h5",Float64,N,N,"w")
#			@test true
#		end
#		@testset "DiskFourTensor" begin
#			DiskFourTensor("/tmp/dtens1.h5",Float64,10,10,10,10,"w")
#			@test true
#		end
#	end
#	@testset "Dot" begin
#		@testset "DiskVector" begin
#			N = 20000
#			x = DiskVector("/tmp/testdv1.h5",Float64,N,"w")
#			y = DiskVector("/tmp/testdv2.h5",Float64,N,"w")
#			blockfill!(x,1.0)
#			blockfill!(y,1.0)
#			@test (abs(dvdot(x,y,5000) - 2E4) < tol) && (abs(dvdot(x,y) - 2E4 < tol))
#			dvwrite!(x,2.0,:)
#			@test abs(dvdot(x,y) - 4E4) < tol
#		end
#		@testset "DiskMatrix" begin
#			N = 1000
#			a = DiskMatrix("/tmp/matrix1.h5",Float64,N,N,"w")
#			b = DiskMatrix("/tmp/matrix2.h5",Float64,N,N,"w")
#			blockfill!(a,1.0)
#			blockfill!(b,2.0)
#			dmdot(a,b,201) - 2E6
#			@test abs(dmdot(a,b,201) - 2E6) < tol
#		end
#		@testset "DiskFourTensor" begin
#			DiskFourTensor("/tmp/dtens1.h5",Float64,10,10,10,10,"w")
#			@test true
#		end
#	end
#end
