module UCCSD

using JuES.Wavefunction
using JuES.Transformation
using TensorOperations
using LinearAlgebra

include("Denominators.jl")

export do_uccsd

function do_uccsd(refWfn::Wfn,maxit; doprint::Bool=false, return_T2::Bool=false)
    nocca = refWfn.nalpha
    noccb = refWfn.nbeta
    nvira = refWfn.nvira
    nvirb = refWfn.nvirb
    N = max(nocca,noccb)
    Δ = nocca - noccb
    uvsr = refWfn.uvsr
    Ca = refWfn.Ca
    Cb = refWfn.Cb
    Cao = refWfn.Cao
    Cbo = refWfn.Cbo
    Cav = refWfn.Cav
    Cbv = refWfn.Cbv
    h = refWfn.hao
    @tensor begin
        ha[p,q] := Ca[m,q]*h[m,n]*Ca[n,p]
    end
    @tensor begin
        hb[p,q] := Cb[m,q]*h[m,n]*Ca[n,p]
    end
    xoxo = permutedims(tei_transform(uvsr,Ca,Ca,Cao,Cao,"xoxo"),[1,3,2,4])
    xxoo = permutedims(tei_transform(uvsr,Ca,Cao,Ca,Cao,"xxoo"),[1,3,2,4])
    XOXO = permutedims(tei_transform(uvsr,Cb,Cb,Cbo,Cbo,"XOXO"),[1,3,2,4])
    fa = form_fa(ha,xoxo,xxoo,XOXO)
    xxoo = nothing
    XXOO = permutedims(tei_transform(uvsr,Cb,Cbo,Cb,Cbo,"XXOO"),[1,3,2,4])
    fb = form_fb(hb,XOXO,XXOO,xoxo)
    xoxo = nothing
    XOXO = nothing
    XXOO = nothing
    
end
function form_Fae(fa,oovv,oOvV,OovV,vovv,vOvV,tia,tiatia,tiatIA,tIAtia,tijab,tiJaB,tIjaB)
end
function form_FAE(fb,OOVV,OoVv,oOVv,VOVV,VoVv,tIA,tIAtIA,tIAtia,tiatIA,tIJAB,tIjAb,tiJAb)
end
function form_Fmi(fa,ooov,oOoV,oovv,oOvV,oOVv,tia,tIA,tiatia,tiatIA,tijab,tiJaB,tiJAb)
end
function form_FMI(fb,OOOV,OoOv,OOVV,OoVv,OovV,tia,tIA,tIAtIA,tIAtia,tIJAB,tIjAb,tIjaB)
end
function form_Fme(fa,oovv,oOvV,tia,tIA)
end
function form_FME(fb,OOVV,OoVv,tia,tIA)
end
function form_Wmnij(oooo,ooov,oovv,tia,tiatia,tijab)
end
function form_WmNiJ(oOoO,oOoV,oOOv,oOvV,oOVv,tia,tIA,tiatIA,tiJaB,tiJAb)
end
end #module
