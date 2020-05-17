using JuES
using Test

#e = JuES.Input.run("ccd_sto3g.dat")
#@test e ≈ -0.07015050066089029
#e = JuES.Input.run("ccsd_sto3g.dat")
e = JuES.Input.run("mccsd_sto3g.dat")
JuES.Input.run("ccsd_sto3g.dat")
JuES.Input.run("mccsd_tz.dat")
@test  isapprox(e,-0.070680102078571;atol=1E-10)

