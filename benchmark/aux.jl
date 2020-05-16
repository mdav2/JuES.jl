using JuES
using JuES.CoupledCluster: RCCSD,RCCD,DFRCCD,mRCCD,mRCCSD
using JLD
wfn = load("/tmp/mywfn.jld","mywfn")
printdo = false
println("RCCSD")
println(mRCCSD.do_rccsd(wfn; doprint=printdo))
