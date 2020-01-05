using Test
using PyCall
using JuES.Wavefunction
using JuES.CISingles
using JuES.CoupledCluster
using JuES.DiskTensors
psi4 = pyimport("psi4")
psi4.core.be_quiet() #turn off output

@testset "All" begin
include("TestWavefunction.jl")
include("TestCISingles.jl")
include("TestDiskTensors.jl")
end
