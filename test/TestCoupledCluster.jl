using JuES

e = JuES.Input.run("ccd_sto3g.dat")
@test e ≈ -0.07015050066089029
e = JuES.Input.run("ccsd_sto3g.dat")
@test e ≈ -0.070680102078571
