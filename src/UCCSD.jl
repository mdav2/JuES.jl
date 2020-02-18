module UCCSD

using JuES.Wavefunction
using JuES.Transformation
using TensorOperations
using LinearAlgebra

include("Denominators.jl")

export do_uccsd

function do_uccsd(refWfn::Wfn; maxit=40, doprint::Bool=false, return_T2::Bool=false)
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
        hb[p,q] := Cb[m,q]*h[m,n]*Cb[n,p]
    end


    oooo = tei_transform(uvsr,Cao,Cao,Cao,Cao,"oooo")
    OOOO = tei_transform(uvsr,Cbo,Cbo,Cbo,Cbo,"OOOO")
    ooOO = tei_transform(uvsr,Cao,Cao,Cbo,Cbo,"ooOO")

    escf = refWfn.vnuc
    rocca = 1:nocca
    roccb = 1:noccb
    for i in rocca
        escf += ha[i,i]
    end
    for i in roccb
        escf += ha[i,i]
    end
    for i in rocca
        for j in rocca
            escf += 0.5*oooo[i,i,j,j]
            escf -= 0.5*oooo[i,j,j,i]
        end
    end
    for i in roccb
        for j in roccb
            escf += 0.5*OOOO[i,i,j,j]
            escf -= 0.5*OOOO[i,j,j,i]
        end
    end
    for i in rocca
        for j in roccb
            escf += ooOO[i,i,j,j]
        end
    end
    if doprint println("@UHF $escf") end

    xxoo = tei_transform(uvsr,Ca,Ca,Cao,Cao,"xxoo")
    xoxo = tei_transform(uvsr,Ca,Cao,Ca,Cao,"xoxo")
    XXOO = tei_transform(uvsr,Cb,Cb,Cbo,Cbo,"XXOO")
    XOXO = tei_transform(uvsr,Cb,Cbo,Cb,Cbo,"XOXO")

    #these dont work, so just creating canonical fock matrix
    #fa = form_fa(ha,xxoo,xoxo,XXOO)
    #fb = form_fb(hb,XXOO,XOXO,xxoo)
    fa = diagm(refWfn.epsa)
    fb = diagm(refWfn.epsb)

    iajb = tei_transform(uvsr,Cao,Cav,Cao,Cav,"iajb")
    oovv = permutedims(iajb,[1,3,2,4]) - permutedims(permutedims(iajb,[1,4,3,2]),[1,3,2,4])
    iajb = nothing
    IAJB = tei_transform(uvsr,Cbo,Cbv,Cbo,Cbv,"IAJB")
    OOVV = permutedims(IAJB,[1,3,2,4]) - permutedims(permutedims(IAJB,[1,4,3,2]),[1,3,2,4])
    oOvV = permutedims(tei_transform(uvsr,Cao,Cav,Cbo,Cbv,"oOvV"),[1,3,2,4])


    Dia   = form_Dia(oovv,fa)
    DIA   = form_Dia(OOVV,fb)
    Dijab = form_Dijab(oovv,fa)
    DIJAB = form_Dijab(OOVV,fb)
    DiJaB = form_DiJaB(oOvV,fa,fb)

    tia = fa[1:nocca,nocca:size(fa)[1]] ./ Dia
    tIA = fb[1:noccb,noccb:size(fb)[1]] ./ DIA
    tijab = oovv ./ Dijab
    tIJAB = OOVV ./ DIJAB
    tiJaB = oOvV ./ DiJaB

    @tensor begin
        emp2_aa[] := (1/4)*tijab[i,j,a,b]*oovv[i,j,a,b]
        emp2_bb[] := (1/4)*tIJAB[i,j,a,b]*OOVV[i,j,a,b]
        emp2_ab[] := tiJaB[i,j,a,b]*oOvV[i,j,a,b]
    end
    emp2 = emp2_aa[] + emp2_bb[] + emp2_ab[]
    if doprint println("@UMP2 $emp2") end
    return emp2

end
function form_fa(ha,xxoo,xoxo,XXOO)
    fa = ha
    R = 1:size(xxoo)[1]
    rocca = 1:size(xxoo)[4]
    roccb = 1:size(XXOO)[4]
    #@tensor begin
    #    fa[p,q] += ha[p,q] + xxoo[p,q,k,k] - xoxo[p,k,q,k] + XXOO[p,q,k,k]
    #end
    for p in R
        for q in R
            tsum = 0
            for k in rocca
                tsum += xxoo[p,q,k,k] - xoxo[p,k,q,k]
            end
            for k in roccb
                tsum += XXOO[p,q,k,k]
            end
            fa[p,q] += tsum
        end
    end
    return fa
end
function form_fb(hb,XXOO,XOXO,xxoo)
    fb = zeros(size(hb))
    @tensor begin
        fb[p,q] = hb[p,q] + XXOO[p,q,k,k] - XOXO[p,k,q,k] + xxoo[p,q,k,k]
    end
    return fb
end
function form_Fae(fa,oovv,oOvV,vovv,vOvV,tia,tiatia,tiatIA,tijab,tiJaB)
    Fae = deepcopy(fa) - diagm(diag(fa))
    @tensor begin
        Fae[a,e] += (-fa[m,e]*tia[m,a]
                    +tia[m,f]*vovv[a,m,e,f]
                    +tIA[m,f]*vOvV[a,m,e,f]
                    -(1/2)*(tijab[m,n,a,f] + (1/2)*(tiatia[m,n,a,f] - tiatia[m,n,f,a]))*oovv[m,n,e,f]
                    -(tiJaB[m,n,a,f] + (1/2)*tiatIA[m,n,a,f])*oOvV[m,n,e,f])
    end
    return Fae
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
