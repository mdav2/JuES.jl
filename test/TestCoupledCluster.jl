using Test
using JuES

JuES.Output.set_print("none")
e_RCCD   = JuES.Input.run("ccd_sto3g.dat")
e_RCCSD  = JuES.Input.run("ccsd_sto3g.dat")
e_mRCCD  = JuES.Input.run("mccd_sto3g.dat")
e_mRCCDs = JuES.Input.run("mccd_sto3g_Float32.dat")

@testset "CoupledCluster" begin
    @testset "Smoke" begin
        @test isapprox(e_RCCD ,   -0.07015050066089029; atol=1E-10)
        @test isapprox(e_mRCCD  , -0.07015050066089029; atol=1E-10)
        @test isapprox(e_mRCCDs , -0.07015049291217820; atol=1E-10)
        @test isapprox(e_RCCSD , -0.070680102078571; atol=1E-10)
        #@test mRCCSD.do_rccsd(JuWfn2) ≈ -0.070680102078571
        #@test RCCSD.do_rccsd(JuWfn2_s) ≈ -0.07068009398829002
    end
end
