module DFRCCD
using JuES.DF
using JuES.Wavefunction
using TensorOperations
using JuES.Transformation
include("Denominators.jl")
export do_df_rccd
"""
    do_df_rccd
"""
function do_df_rccd(refWfn::Wfn; maxit=40, doprint=false, return_T2=false, dfbname="default")
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    epsa = refWfn.epsa
    T = eltype(refWfn.uvsr)
    bov,bvo,boo,bvv = make_df_rccd_integrals(refWfn::Wfn; dfbname=dfbname)
    T2 = zeros(T, nocc, nocc, nvir, nvir)
    Dijab = form_Dijab(T2, epsa)
    T2_init!(T2, bov, Dijab)
    ecc = ccenergy(T2, bov)
    if doprint println("@DF-RMP2 $ecc") end
    Fae = form_Fae(T2, bov)
    return ecc
end

function make_df_rccd_integrals(refWfn::Wfn; dfbname="default")
    bμν = DF.make_bμν(refWfn; dfbname=dfbname)
    Cao = refWfn.Cao
    Cav = refWfn.Cav
    @tensoropt begin
        bov[i,a,Q] := Cao[μ,i]*bμν[μ,ν,Q]*Cav[ν,a]
        bvo[a,i,Q] := Cav[μ,a]*bμν[μ,ν,Q]*Cao[ν,i]
        boo[i,j,Q] := Cao[μ,i]*bμν[μ,ν,Q]*Cao[ν,j]
        bvv[a,b,Q] := Cav[μ,a]*bμν[μ,ν,Q]*Cav[ν,b]
    end
    return bov,bvo,boo,bvv
end
function ccenergy(tiJaB, bov)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    for i in rocc
        for j in rocc
            bi = bov[i,:,:]
            bj = bov[j,:,:]
            @tensor begin
                bAB[a,b] := bi[a,Q]*bj[b,Q]
            end
            for a in rvir
                for b in rvir
                    ecc += bAB[a,b] * 2 * tiJaB[i,j,a,b]
                    ecc -= bAB[a,b] * tiJaB[j,i,a,b]
                end
            end
        end
    end
    return ecc
end
function T2_init!(tiJaB,bov,Dijab)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    nao  = size(bov,   3)
    rocc = 1:nocc
    rvir = 1:nvir
    #@tensor begin
    #    tiJaB[i,j,a,b] = bov[i,a,Q]*bov[j,b,Q]#/Dijab[i,j,a,b]
    #end
    tiJaB .= tiJaB ./ Dijab
    for i in rocc
        for j in rocc
            bi = bov[i,:,:]
            bj = bov[j,:,:]
            @tensor begin
                bAB[a,b] := bi[a,Q]*bj[b,Q]
            end
            for a in rvir
                for b in rvir
                    tiJaB[i,j,a,b] = bAB[a,b]/Dijab[i,j,a,b]
                end
            end
        end
    end
end
function form_Fae(tiJaB, bov)
    dt = eltype(bov)
    nvir = size(tiJaB, 4)
    Fae = zeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, bov)
    return Fae
end
function form_Fae!(Fae,tiJaB, bov)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = 1:nocc
    rvir = 1:nvir
    @tensor begin
        iajb[i,a,j,b] := bov[i,a,Q]*bov[j,b,Q]
    end
    #@tensoropt begin
    #    Fae[a,e] = -1*bov[m,e,Q]*bov[n,f,Q]*(2*tiJaB[m,n,a,f] - tiJaB[n,m,a,f])
    #end
    @tensoropt begin
        Fae[a,e] = -1*iajb[m,e,n,f]*(2*tiJaB[m,n,a,f] - tiJaB[n,m,a,f])
    end
    return Fae
end
function form_Fmi(tiJaB, bov)
    dt = eltype(bov)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Fmi = zeros(dt, nocc, nocc)
    form_Fmi!(Fmi, tiJaB, bov)
    return Fmi
end
function form_Fmi!()
end
end #module
