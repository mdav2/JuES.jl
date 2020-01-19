module ROCCD

using JuES.Wavefunction
using JuES.Transformation
using TensorOperations
include("Denominators.jl")

export do_roccd

function do_roccd(refWfn::Wfn,maxit; doprint::Bool=false, return_T2::Bool=false)
    nocca = refWfn.nalpha
    noccb = refWfn.nbeta
    nvira = refWfn.nvira
    nvirb = refWfn.nvirb
    uvsr = refWfn.uvsr
    Ca = refWfn.Ca
    Cb = refWfn.Cb
    Cao = refWfn.Cao
    Cbo = refWfn.Cbo
    Cav = refWfn.Cav
    Cbv = refWfn.Cbv
    h = refWfn.hao
    @tensor begin
        ha[p,q] := Ca[p,μ]*h[μ,ν]*Ca[ν,q]
    end
    @tensor begin
        hb[p,q] := Cb[p,μ]*h[μ,ν]*Ca[ν,q]
    end
    xoxo = permutedims(tei_transform(uvsr,Ca,Ca,Cao,Cao,"xoxo")[1,3,2,4])
    xxoo = permutedims(tei_transform(uvsr,Ca,Cao,Ca,Cao,"xxoo")[1,3,2,4])
    XOXO = permutedims(tei_transform(uvsr,Cb,Cb,Cbo,Cbo,"XOXO")[1,3,2,4])
    fa = form(fa,ha,xoxo,xxoo,XOXO)
    xxoo = nothing
    XXOO = permutedims(tei_transform(uvsr,Cb,Cbo,Cb,Cbo,"XXOO")[1,3,2,4])
    fb = form(fb,hb,XOXO,XXOO,xoxo)
    xoxo = nothing
    XOXO = nothing
    XXOO = nothing

    oovv = permutedims(tei_transform(uvsr,Cao,Cav,Cao,Cav,"oovv"),[1,3,2,4])
    oOvV = permutedims(tei_transform(uvsr,Cao,Cav,Cbo,Cbv,"oOvV"),[1,3,2,4])
    OOVV = permutedims(tei_transform(uvsr,Cbo,Cbv,Cbo,Cbv,"OOVV"),[1,3,2,4])
    vVvV = permutedims(tei_transform(uvsr,Cav,Cav,Cbv,Cbv,"vVvV"),[1,3,2,4])
    oVvO = permutedims(tei_transform(uvsr,Cao,Cav,Cbv,Cbo,"oVvO"),[1,3,2,4])
    vOvO = permutedims(tei_transform(uvsr,Cav,Cav,Cbo,Cbo,"vOvO"),[1,3,2,4])
    oVoV = permutedims(tei_transform(uvsr,Cao,Cao,Cbv,Cbv,"oVoV"),[1,3,2,4])
    vOoV = permutedims(tei_transform(uvsr,Cav,Cao,Cbo,Cbv,"vOoV"),[1,3,2,4])

    Dijab = form_Dijab(oovv,fa)
    DiJaB = form_DiJaB(oOvV,fa,fb)
    DIJAB = form_Dijab(OOVV,fb)

    tijab = zeros(size(oovv))
    tiJaB = zeros(size(oOvV))
    tIJAB = zeros(size(OOVV))

    tijab .= oovv ./ Dijab
    tiJaB .= oOvV ./ DiJaB
    tIJAB .= OOVV ./ DIJAB

    Fae = zeros(nvira,nvira)
    FAE = zeros(nvirb,nvirb)
    Fmi = zeros(nocca,nocca)
    FMI = zeros(noccb,noccb)
    Fme = fa[1:nocca,nocca+1:]
    FME = fb[1:noccb,noccb+1:]
    Wmnij = zeros(nocca,nocca,nocca,nocca)
    WmNiJ = zeros(nocca,noccb,nocca,noccb)
    WMNIJ = zeros(noccb,noccb,noccb,noccb)
    Wabef = zeros(nvira,nvira,nvira,nvira)
    WaBeF = zeros(nvira,nvirb,nvira,nvirb)
    WABEF = zeros(nvirb,nvirb,nvirb,nvirb)
    Wmbej = zeros(nocca,nvira,nvira,nocca)
    WmBeJ = zeros(nocca,nvirb,nvira,noccb)
    WmBEj = zeros(nocca,nvirb,nvirb,nocca)
    WMBEJ = zeros(noccb,nvirb,nvirb,noccb)
    WMbEj = zeros(noccb,nvira,nvirb,nocca)
    WMbeJ = zeros(noccb,nvira,nvira,noccb)


    form_Fae!(Fae,fa[nocca+1:,nocca+1:],oOvV,tijab,tiJaB)
    form_Fae!(FAE,fb[noccb+1:,noccb+1:],oOvV,tIJAB,tiJaB)
    form_Fmi!(Fmi,fa[1:nocca,1:nocca],oOvV,tijab,tiJaB)
    form_Fmi!(FMI,fb[1:noccb,1:noccb],oOvV,tIJAB,tiJaB)
    form_Wmnij!(Wmnij,oOvV,tijab)
    form_WmNiJ!(WmNiJ,oOvV,tiJaB)
    form_Wmnij!(WMNIJ,oOvV,tIJAB)
    form_Wabef!(Wabef,vVvV,oOvV,tijab)
    form_WaBeF!(WaBeF,vVvV,oOvV,tiJaB)
    form_WABEF!(WABEF,vVvV,oOvV,tIJAB)
    form_Wmbej!(Wmbej,oVvO,vOvO,oOvV,tijab,tiJaB)
    form_WmBeJ!(WmBeJ,oVvO,oOvV,tIJAB,tiJaB)
    form_WmBEj!(WmBEj,oVoV,oOvV,tiJaB)
    form_WMBEJ!(WMBEJ,oVvO,vOvO,oOvV,tIJAB,tiJaB)
    form_WMbEj!(WMbEj,vOoV,oOvV,tijab,tiJaB)
    form_WMbeJ!(WMbeJ,vOvO,oOvV,tiJaB)

end
function form_tijab(oovv,Fae,Fmi,Wmnij,Wabef,Wmbej,WMbEj,tijab,tiJaB,Dijab)
    @tensoropt begin
        _tijab[i,j,a,b] := (oovv[i,j,a,b]
                            + tijab[i,j,a,e]*Fae[b,e]
                            - tijab[i,j,b,e]*Fae[a,e]
                            - tijab[i,m,a,b]*Fmi[m,j]
                            + tijab[j,m,a,b]*Fmi[m,i]
                            + (1/2)*tijab[m,n,a,b]*Wmnij[m,n,i,j]
                            + (1/2)*tijab[i,j,e,f]*Wabef[a,b,e,f]
                            + tijab[i,m,a,e]*Wmbej[m,b,e,j]
                            + tiJaB[i,m,a,e]*WMbEj[m,b,e,j]
                            - tijab[i,m,b,e]*Wmbej[m,a,e,j]
                            - tiJaB[i,m,b,e]*WMbEj[m,a,e,j]
                            - tijab[j,m,a,e]*Wmbej[m,b,e,i]
                            - tiJaB[j,m,a,e]*WMbEj[m,b,e,i]
                            + tijab[j,m,b,e]*Wmbej[m,a,e,i]
                            + tiJaB[j,m,b,e]*WMbEj[m,a,e,i])
    end
    return _tijab
end
#function form_tiJaB(oOvV,Fae,FAE,Fmi,FMI,WmN
function form_fa(ha,xoxo,xxoo,XOXO)
    fa = zeros(size(ha))
    @tensor begin
        fa[p,q] = ha[p,q] + xoxo[p,k,q,k] - xxoo[p,q,k,k] + XOXO[p,K,q,K]
    end
end
function form_fb(hb,XOXO,XXOO,xoxo)
    fb = zeros(size(hb))
    @tensor begin
        fb[p,q] = hb[p,q] + XOXO[p,k,q,k] - XXOO[p,q,k,k] + xoxo[p,K,q,K]
    end
end
function form_Fae!(Fae,fa,oOvV,tijab,tiJaB)
    @tensoropt begin
        Fae[a,e] := -0.5*(tijab[m,n,a,f]*(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
                    - tiJaB[m,n,a,f]*oOvV[m,n,e,f]
    end
    Fae += fa - diagm(diag(fa))
    return Fae
end
function form_Fmi!(Fmi,fa,oOvV,tijab,tiJaB)
    @tensoropt begin
        Fmi[m,i] := 0.5*(tijab[i,n,e,f]*(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
                    + tiJaB[i,n,e,f]*oOvV[m,n,e,f]
    end
    Fmi += fa - diagm(diag(fa))
end
function form_Wmnij!(Wmnij,oOvV,tijab)
    @tensoropt begin
        Wmnij[m,n,i,j] = (oOvV[m,n,i,j] - oOvV[n,m,i,j]
                          + (1/4)*(tijab[i,j,e,f]
                                   *(oOvV[m,n,e,f] - oOvV[m,n,f,e])))
    end
    return Wmnij
end
function form_WmNiJ!(WmNiJ,oOvV,tiJaB)
    @tensoropt begin
        WmNiJ[m,n,i,j] = (oOvV[m,n,i,j] + (1/4)*tiJaB[i,j,e,f]*oOvV[m,n,e,f]
                          + (1/4)*tiJab[i,j,f,e]*oOvV[m,n,f,e])
    end
    return WmNiJ
end
function form_Wabef!(Wabef,vVvV,oOvV,tijab)
    @tensoropt begin
        Wabef[a,b,e,f] = (vVvV[a,b,e,f] - vVvV[b,a,e,f] 
                          + (1/4)*tijab[m,n,a,b]
                          *(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
    end
end
function form_WaBeF!(WaBeF,vVvV,oOvV,tiJaB)
    @tensoropt begin
        WaBeF[a,b,e,f] = vVvV[a,b,e,f] + (1/2)*tiJaB[m,n,a,b]*oOvV[m,n,e,f]
    end
end
function form_Wmbej!(Wmbej,oVvO,vOvO,oOvV,tijab,tiJaB)
    @tensoropt begin
        Wmbej[m,b,ej] = (oVvO[m,b,e,j] - vOvO[b,m,e,j] 
                         - (1/2)*tijab[j,n,f,b]
                         *(oOvV[m,n,e,f] - oOvV[n,m,e,f])
                         + (1/2)*tiJaB[j,n,b,f]*oOvV[m,n,e,f])
    end
    return Wmbej
end
function form_WmBeJ!(WmBeJ,oVvO,oOvV,tIJAB,tiJaB)
    @tensoropt begin
        WmBeJ[m,b,e,j] = (oVvO[m,b,e,j] 
                          - (1/2)*tIJAB[j,n,f,b]*oOvV[m,n,e,f]
                          + (1/2)*tiJaB[n,j,f,b]
                            *(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
    end
    return WmBeJ
end
function form_WmBEj!(WmBEj,oVoV,oOvV,tiJaB)
    @tensoropt begin
        WmBEj[m,b,e,j] = (-1*oVoV[m,b,j,e]
                          + (1/2)*tiJaB[j,n,f,b]*oOvV[m,n,f,e])
    end
    return WmBEj
end
function form_WMBEJ!(WMBEJ,oVvO,vOvO,oOvV,tIJAB,tiJaB)
    @tensoropt begin
        WMBEJ[m,b,e,j] = (oVvO[m,b,e,j] - vOvO[b,m,e,j]
                          - (1/2)*tIJAB[j,n,f,b]
                            *(oOvV[m,n,e,f] - oOvV[n,m,e,f])
                          + (1/2)*tiJaB*oOvV[n,m,f,e])
    end
    return WMBEJ
end
function form_WMbEj!(WMbEj,vOoV,oOvV,tijab,tiJaB)
    @tensoropt begin
        WMbEj[m,b,e,j] = (vOoV[b,m,j,e]
                          - (1/2)*tijab[j,n,f,b]*oOvV[n,m,f,e]
                          + (1/2)*tiJaB[j,n,b,f]
                            *(oOvV[m,n,e,f] - oOvV[n,m,e,f]))
    end
    return WMbEj
end
function form_WMbeJ(WMbeJ,vOvO,oOvV,tiJaB)
    @tensoropt begin
        WMbeJ[m,b,e,j] = (-1*vOvO[b,m,e,j]
                          + (1/2)*tiJaB[n,j,b,f]*oOvV[n,m,e,f])
    end
    return WMbeJ
end
