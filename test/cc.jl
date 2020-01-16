using Test
using JuES
using Base.Threads


psi4.core.be_quiet() #turn off output
tol = 1E-14
mol2 = psi4.geometry("""
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("basis" => "cc-pvdz", "scf_type" => "pk", "d_convergence" => 14))
e, wfn2 = psi4.energy("hf/sto-3g", mol = mol2, return_wfn = true)
println("echo")
Threads.@spawn println("echo")
JuWfn2 = JuES.Wavefunction.Wfn(wfn2)
JuWfn3 = JuES.Wavefunction.Wfn(wfn2, Float64, true, true)
@testset "CoupledCluster" begin
    @testset "Smoke" begin
#        Threads.@spawn println("echo")
        @test JuES.CoupledCluster.do_rccd(JuWfn2, 40, doprint=false) ≈ -0.07015050066089029
        @test JuES.CoupledCluster.do_rccd(JuWfn3, 40, doprint=false) ≈ -0.07015050066089029
    end
end
