using JuES
using JuES.CoupledCluster.mRCCSD

JuES.Psi4.psi4.core.be_quiet() #turn off output
JuES.Psi4.println(Threads.nthreads())
JuES.Psi4.psi4.set_num_threads(6)
JuES.Psi4.psi4.set_options(Dict("D_CONVERGENCE" => 14,
                      "E_CONVERGENCE" => 12,
					  "scf_type" => "pk"))
mol2 = JuES.Psi4.psi4.geometry("""
                     O
                     H 1 1.1
                     H 1 1.1 2 104.0
					 symmetry c1
					 """)
e2,wfn2 = JuES.Psi4.psi4.energy("hf/cc-pvtz",mol=mol2,return_wfn=true)
JuWfn2 = Wfn(wfn2)

