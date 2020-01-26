"""
Test suite for JuES. 
To run tests:
	shell> julia
	julia> import Pkg
	julia> Pkg.test("JuES")

to add a test:
1) determine a reliable and simple test case!
2) write a self contained file called Test<Descriptive/ModuleName>.jl
	tests should be organized like
	shell> cat Test<Descriptive/ModuleName>.jl
	<imports>
	<setup>
	@testset "Descriptive/ModuleName" begin
		@testset "subset1" begin
			@test true
		end
	end

2) add a line in testset All (this file) calling your tests
	@testset "All" begin
		<old tests>
		include("<Descriptive/ModuleName>.jl")
	end
"""

using Test
using PyCall
using JuES.Wavefunction
using JuES.CISingles
using JuES.CoupledCluster
using JuES.DiskTensors
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output

@testset "JuES" begin
    include("TestWavefunction.jl")
    include("TestCISingles.jl")
    include("TestDiskTensors.jl")
    include("TestCoupledCluster.jl")
    include("TestTransformation.jl")
    include("TestMollerPlesset.jl")
end
