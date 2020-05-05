module CISingles

using JuES.Wavefunction
using JuES.Transformation
using JuES

function do_RCIS(refWfn::Wfn; dotriplets=false)
    Cao = refWfn.Cao
    Cav = refWfn.Cav
    vovo = tei_transform(refWfn.uvsr,Cav,Cao,Cav,Cao,"vovo")
    vvoo = tei_transform(refWfn.uvsr,Cav,Cav,Cao,Cao,"vvoo")
end

end#module
