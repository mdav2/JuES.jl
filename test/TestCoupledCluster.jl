using Test
using PyCall
using JuES
using JuES.Wavefunction
using JuES.CoupledCluster: RCCSD, RCCD, ROCCD


psi4.core.be_quiet() #turn off output
psi4.set_num_threads(6)
tol = 1E-14
mol2 = psi4.geometry("""
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("basis" => "sto-3g", "scf_type" => "pk", "d_convergence" => 14))
e, wfn2 = psi4.energy("hf/sto-3g", mol = mol2, return_wfn = true)
JuWfn2 = Wfn(wfn2)
mol3 = psi4.geometry("""
      1 2
      O
      H 1 1.1
      H 1 1.1 2 104.0
      symmetry c1
      """)
psi4.set_options(Dict("reference" => "rohf"))
e, wfn3 = psi4.energy("hf/sto-3g",mol=mol3,return_wfn=true)
JuWfn3 = Wfn(wfn3)
#JuWfn3 = Wfn(wfn2, Float64, true, true) #disk based CCD is currently NOT working
@testset "CoupledCluster" begin
    @testset "Smoke" begin
        @test RCCD.do_rccd(JuWfn2, 40, doprint=false) ≈ -0.07015050066089029
        #@test do_rccd(JuWfn3, 40, doprint=false) ≈ -0.07015050066089029
        @test RCCSD.do_rccsd(JuWfn2, 40, doprint=false) ≈ -0.070680102078571
        ROCCD.do_roccd(JuWfn3, 40, doprint=true)
    end
end
