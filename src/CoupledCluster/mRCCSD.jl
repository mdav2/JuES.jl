"""
    JuES.CoupledCluster.RCCSD

## Functions

    JuES.CoupledCluster.RCCSD.do_rccsd
"""
module mRCCSD
using JuES.Wavefunction
using JuES.Transformation
using JuES.IntegralTransformation
using JuES.Output
using Printf
using LinearAlgebra
using JuES
using TensorOperations
include("Denominators.jl")

"""
    JuES.CoupledCluster.RCCSD.do_rccsd(refWfn::Wfn; kwargs...)

## Arguments
    refWfn::Wfn
reference wavefunction to compute RCCSD on.

    kwargs::Any
keyword arguments for various options.

## Kwargs
    :doprint::Bool  print or no? deprecated
    :maxit::Int     max number of iterations for RCCSD equations
    :return_T::Bool do return the converged cluster amplitudes?
    :diis::Bool     do DIIS extrapolation? currently dummy  
"""
function do_rccsd(refWfn::Wfn; kwargs...)
    for arg in keys(JuES.CoupledCluster.defaults)
        if arg in keys(kwargs)
            @eval $arg = $(kwargs[arg])
        else
            @eval $arg = $(JuES.CoupledCluster.defaults[arg])
        end
    end
    JuES.CoupledCluster.print_header()
    @output "*   executing CCSD\n"
    maxit = 40
    return_T = false
    diis = false
    set_zero_subnormals(true)
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    uvsr = refWfn.ao_eri
    Cav = refWfn.Cav
    Cao = refWfn.Cao
    dtt = eltype(uvsr)

    @output "    Forming MO basis integrals ... "
    t = @elapsed begin
        oooo, ooov, oovv, ovov, ovvv, vvvv = make_rccsd_integrals(refWfn)
    end
    @output "done in {:>5.2f}s\n" t
    T_total = @elapsed begin
        epsa = refWfn.epsa
        T2 = zeros(dtt, nocc, nocc, nvir, nvir)
        T1 = zeros(dtt,nocc,nvir)
        Fae = zeros(nvir,nvir)
        Fmi = zeros(nocc,nocc)
        Fme = zeros(nocc,nvir)
        Wabef = zeros(nvir,nvir,nvir,nvir)
        __Wabef = zeros(nvir,nvir*nvir*nvir)
        __tijab = zeros(nocc,nvir,nocc,nvir)
        __oovv = zeros(nocc,nvir,nocc,nvir)
        Wmnij = zeros(nocc,nocc,nocc,nocc)
        WmBeJ = zeros(nocc,nvir,nvir,nocc)
        WmBEj = zeros(nocc,nvir,nvir,nocc)
        Dia = [epsa[i] - epsa[a] for i=1:nocc,a=nocc+1:nocc+nvir]
        Dijab = [epsa[i] + epsa[j] - epsa[a] - epsa[b] for i=1:nocc,j=1:nocc,a=nocc+1:nocc+nvir,b=nocc+1:nocc+nvir]
        T2 .= oovv ./ Dijab
        o = nocc
        v = nvir
        _T1 = reshape(T1,o*v)
        tiatia = zeros(dtt,o*v,o*v)
        BLAS.gemm!('N','T',1.0,_T1,_T1,0.0,tiatia)
        tiatia = permutedims(reshape(tiatia,o,v,o,v),(1,3,2,4))


        @output "    @MP2 {:>20.17f}\n" ccenergy(oovv,T1,tiatia,T2)
        @output "    Starting CC iterations\n\n"
        @output "{:>7s} {:>20}\n" "Iter." "E[CCSD]"
        @output repeat("~",80)*"\n"
        for i in 1:maxit
            _T1,_T2 = cciter(Fae,Fmi,Fme,
                   Wabef,Wmnij,WmBeJ,WmBEj,
                   oooo,ooov,oovv,ovov,ovvv,vvvv,
                   T1,tiatia,T2,Dia,Dijab,__Wabef,__tijab,__oovv)
            T1 = _T1
            T2 = _T2
            @tensor begin
                tiatia[m,n,a,f] = T1[m,a]*T1[n,f]
            end
            output("    {:<3d} {:>20.17f}\n",i,ccenergy(oovv,T1,tiatia,T2))

        end
    end
    output(repeat("~",80)*"\n")
    @output "    CC iterations completed in {:>5.2f}\n" T_total
    @output "    CC final energy:\n"
    @output "        @CCSD {:>20.17f}\n" ccenergy(oovv,T1,tiatia,T2)
    if return_T
        return ccenergy(oovv,T1,tiatia,T2),T2
    else
        return ccenergy(oovv,T1,tiatia,T2)
    end
end
function make_rccsd_integrals(wfn)
    oooo = get_eri(wfn,"oooo")
    ooov = get_eri(wfn,"ooov")
    oovv = get_eri(wfn,"oovv")
    ovov = get_eri(wfn,"ovov")
    ovvv = get_eri(wfn,"ovvv")
    vvvv = get_eri(wfn,"vvvv")
    return oooo, ooov, oovv, ovov, ovvv, vvvv
end
function ccenergy(ijab,tia,tiatia,tijab)
    @tensor begin
        ecc[] := ijab[i,j,a,b]*(2*tijab[i,j,a,b] + 2*tiatia[i,j,a,b] 
                              - tijab[j,i,a,b] - tiatia[j,i,a,b])
    end
    return ecc[]
end
function cciter(Fae,Fmi,Fme,
                Wabef,Wmnij,WmBeJ,WmBEj,
                oooo,ooov,oovv,ovov,ovvv,vvvv,
                tia,tiatia,tijab,Dia,Dijab,__Wabef,__tijab,__oovv)
    form_Fae!(Fae,ovvv,oovv,tia,tiatia,tijab)
    form_Fmi!(Fmi,ooov,oovv,tia,tiatia,tijab)
    form_Fme!(Fme,oovv,tia)
    form_Wabef!(__Wabef,Wabef,vvvv,ovvv,oovv,tia,tiatia,tijab)
    form_Wmnij!(Wmnij,oooo,ooov,oovv,tia,tiatia,tijab)
    form_WmBeJ!(WmBeJ,ovvv,ooov,oovv,tia,tiatia,tijab,__tijab,__oovv)
    form_WmBEj!(WmBEj,ovov,ovvv,ooov,oovv,tia,tiatia,tijab)
    _tia = form_T1(Fae,Fmi,Fme,oovv,ovov,ooov,ovvv,tia,tijab,Dia)
    _tijab = form_T2(Fae,Fmi,Fme,Wabef,Wmnij,WmBeJ,WmBEj,oovv,
                    ovov,ovvv,ovvv,ooov,tia,tiatia,tijab,Dijab)
    #tia = nothing
    #tijab = nothing
    return _tia,_tijab
end
                
function form_Fae!(Fae,ovvv,oovv,tia,tiatia,tijab)
    o = size(tijab,1)
    v = size(tijab,3)
    @tensor begin
        _ovvv[1,3,2,4] := 2*ovvv[1,2,3,4]
        _ovvv[1,4,3,2] -= ovvv[1,2,3,4]
        _tijab[3,1,2,4] := tijab[1,2,3,4] + 0.5*tiatia[1,2,3,4]
        _oovv[1,2,4,3] := 2*oovv[1,2,3,4]
        _oovv[2,1,4,3] -= oovv[1,2,3,4]
    end
    _ovvv = reshape(_ovvv,o*v,v*v)
    #_ovvv = 2*reshape(permutedims(ovvv,(1,3,2,4)),o*v,v*v) - reshape(permutedims(ovvv,(1,4,3,2)),o*v,v*v)
    _tia = reshape(tia,1,o*v)
    _Fae = reshape(Fae,1,v*v)
    BLAS.gemm!('N','N',1.0,_tia,_ovvv,0.0,_Fae)
    Fae .= reshape(_Fae,v,v)
    #_tijab = reshape(permutedims(tijab,(3,1,2,4)),v,o^2*v) + 0.5*reshape(permutedims(tiatia,(3,1,2,4)),v,o^2*v)
    #@tensor 
    _tijab = reshape(_tijab,v,o^2*v)
    #_oovv = 2*reshape(permutedims(oovv,(1,2,4,3)),o^2*v,v) - reshape(permutedims(oovv,(2,1,4,3)),o^2*v,v)
    #@tensor begin
    #    _oovv[1,2,4,3] := 2*oovv[1,2,3,4]
    #    _oovv[2,1,4,3] -= oovv[1,2,3,4]
    #end
    _oovv = reshape(_oovv,o^2*v,v)
    BLAS.gemm!('N','N',-1.0,_tijab,_oovv,1.0,Fae)
    return Fae
end
function form_Fmi!(Fmi,ooov,oovv,tia,tiatia,tijab)
    o = size(tijab,1)
    v = size(tijab,3)
    @tensor begin
        _ooov[2,4,1,3] := 2*ooov[1,2,3,4]
        _ooov[1,4,2,3] -= ooov[1,2,3,4]
        _oovv[1,2,3,4] := 2*oovv[1,2,3,4]
        _oovv[1,2,4,3] -= oovv[1,2,3,4]
        _tijab[2,3,4,1] := tijab[1,2,3,4]
        _tijab[2,3,4,1] += 0.5*tiatia[1,2,3,4]
    end
    _ooov = reshape(_ooov,o*v,o^2)
    _tia = reshape(tia,1,o*v)
    _Fmi = reshape(Fmi,1,o^2)
    BLAS.gemm!('N','N',1.0,_tia,_ooov,0.0,_Fmi)
    Fmi .= reshape(_Fmi,o,o)
    @tensor begin
    end
    #_oovv = 2*reshape(oovv,o,v^2*o) - reshape(permutedims(oovv,(1,2,4,3)),o,v^2*o)
    _oovv = reshape(_oovv,o,v^2*o)
    #_tijab = reshape(permutedims(tijab,(2,3,4,1)),v^2*o,o) + 0.5*reshape(permutedims(tiatia,(2,3,4,1)),v^2*o,o)
    _tijab = reshape(_tijab,v^2*o,o) 
    BLAS.gemm!('N','N',1.0,_oovv,_tijab,1.0,Fmi)
    return Fmi
end
function form_Fme!(Fme,oovv,tia)
    o = size(oovv,1)
    v = size(oovv,3)
    _oovv = 2*reshape(permutedims(oovv,(1,3,2,4)),o*v,o*v) - reshape(permutedims(oovv,(2,3,1,4)),o*v,o*v)
    _tia = reshape(tia,o*v)
    _Fme = reshape(Fme,o*v)
    BLAS.gemm!('N','N',1.0,_oovv,_tia,0.0,_Fme)
    Fme = reshape(Fme,o,v)
    return Fme
end
function form_Wmnij!(Wmnij,oooo,ooov,oovv,tia,tiatia,tijab)
    o = size(oovv,1)
    v = size(oovv,3)
    Wmnij .= oooo
    _Wmnij = reshape(Wmnij,o^3,o)
    _ooov = reshape(ooov,o^3,v)
    BLAS.gemm!('N','T',1.0,_ooov,tia,1.0,_Wmnij)
    Wmnij .= reshape(_Wmnij,o,o,o,o)

    @tensor _Wmnij[2,1,4,3] := Wmnij[1,2,3,4]
    #_Wmnij = reshape(permutedims(Wmnij,(2,1,4,3)),o^3,o)
    _Wmnij = reshape(_Wmnij,o^3,o)
    BLAS.gemm!('N','T',1.0,_ooov,tia,1.0,_Wmnij)
    _Wmnij = reshape(_Wmnij,o,o,o,o)
    @tensor Wmnij[2,1,4,3] = _Wmnij[1,2,3,4]
    #Wmnij .= permutedims(reshape(_Wmnij,o,o,o,o),(2,1,4,3))

    _tijab = reshape(tijab,o^2,v^2) + reshape(tiatia,o^2,v^2)
    _oovv = reshape(oovv,o^2,v^2)
    _Wmnij = reshape(Wmnij,o^2,o^2)
    BLAS.gemm!('N','T',0.5,_oovv,_tijab,1.0,_Wmnij)
    Wmnij = reshape(Wmnij,o,o,o,o)
    return Wmnij
end
function form_Wabef!(__Wabef,Wabef,vvvv,ovvv,oovv,tia,tiatia,tijab)
    o = size(oovv,1)
    v = size(oovv,3)
    Wabef .= vvvv

    
    _ovvv = reshape(ovvv,o,v^3)
    BLAS.gemm!('T','N',-1.0,tia,_ovvv,0.0,__Wabef)
    _Wabef = reshape(__Wabef,v,v,v,v)
    
    @tensor Wabef[2,1,4,3] += _Wabef[1,2,3,4]
    _Wabef = reshape(Wabef,v,v^3)
    _ovvv = reshape(ovvv,o,v^3)
    BLAS.gemm!('T','N',-1.0,tia,_ovvv,1.0,_Wabef)
    Wabef = reshape(_Wabef,v,v,v,v)
    
    _tijab = tijab + tiatia
    _tijab = reshape(_tijab,o^2,v^2)
    _oovv = reshape(oovv,o^2,v^2)
    _Wabef = reshape(Wabef,v^2,v^2)
    BLAS.gemm!('T','N',0.5,_tijab,_oovv,1.0,_Wabef)
    Wabef = reshape(_Wabef,v,v,v,v)

    #@tensoropt begin
    #    Wabef[a,b,e,f] = (vvvv[a,b,e,f] 
    #                      - tia[m,b]*ovvv[m,a,f,e]
    #                      - tia[m,a]*ovvv[m,b,e,f]
    #                      + 0.5*(tijab[m,n,a,b] + tiatia[m,n,a,b])*oovv[m,n,e,f])
    #end
    return Wabef
end
function form_WmBeJ!(WmBeJ,ovvv,ooov,oovv,tia,tiatia,tijab,__tijab,__oovv)
    o = size(oovv,1)
    v = size(oovv,3)
    @tensor WmBeJ[m,b,e,j] = oovv[m,j,e,b]
    _ovvv = reshape(ovvv,o*v^2,v)
    _WmBeJ = reshape(WmBeJ,o*v^2,o)
    BLAS.gemm!('N','T',1.0,_ovvv,tia,1.0,_WmBeJ)
    WmBeJ = reshape(WmBeJ,o,v,v,o)
    @tensor _WmBeJ[b,m,j,e] := WmBeJ[m,b,e,j]
    _WmBeJ = reshape(_WmBeJ,v,o*o*v)
    _ooov = reshape(ooov,o,o*o*v)
    BLAS.gemm!('T','N',-1.0,tia,_ooov,1.0,_WmBeJ)
    _WmBeJ = reshape(_WmBeJ,v,o,o,v)
    @tensor WmBeJ[m,b,e,j] = _WmBeJ[b,m,j,e]
    @tensor begin
        __tijab[n,f,j,b] = tijab[j,n,f,b] + 2*tiatia[j,n,f,b]
        __oovv[m,e,n,f] = oovv[m,n,e,f]
        _WmBeJ[m,e,j,b] := WmBeJ[m,b,e,j]
    end
    _tijab = reshape(__tijab,o*v,o*v)
    _oovv = reshape(__oovv,o*v,o*v)
    _WmBeJ = reshape(_WmBeJ,o*v,o*v)
    BLAS.gemm!('N','N',-0.5,_oovv,_tijab,1.0,_WmBeJ)
    _WmBeJ = reshape(_WmBeJ,o,v,o,v)
    @tensor WmBeJ[m,b,e,j] = _WmBeJ[m,e,j,b]
    @tensor begin
        WmBeJ[m,b,e,j] += (#oovv[m,j,e,b] 
                          #+ tia[j,f]*ovvv[m,b,e,f]
                          #- tia[n,b]*ooov[n,m,j,e]
                          #- 0.5*(tijab[j,n,f,b] + 2*tiatia[j,n,f,b])*oovv[m,n,e,f]
                          + 0.5*tijab[n,j,f,b]*(2*oovv[m,n,e,f] - oovv[n,m,e,f]))
    end
    return WmBeJ
end
function form_WmBEj!(WmBEj,ovov,ovvv,ooov,oovv,tia,tiatia,tijab)
    @tensor begin
        WmBEj[m,b,e,j] = (-1*ovov[m,b,j,e] - tia[j,f]*ovvv[m,b,f,e]
                          + tia[n,b]*ooov[m,n,j,e]
                          + (0.5*tijab[j,n,f,b] + tiatia[j,n,f,b])*oovv[n,m,e,f])
    end
    return WmBEj
end
function form_T1(Fae,Fmi,Fme,oovv,ovov,ooov,ovvv,tia,tijab,Dia)
    @tensor begin
        _tia[i,a] := (tia[i,e]*Fae[a,e] - tia[m,a]*Fmi[m,i]
                     + Fme[m,e]*(2*tijab[i,m,a,e] - tijab[m,i,a,e])
                     + tia[m,e]*(2*oovv[i,m,a,e] - ovov[m,a,i,e])
                     - tijab[m,n,a,e]*(2*ooov[m,n,i,e] - ooov[n,m,i,e])
                     + tijab[i,m,e,f]*(2*ovvv[m,a,f,e] - ovvv[m,a,e,f]))
    end
    _tia .= _tia ./ Dia
    return _tia
end
function form_T2(Fae,Fmi,Fme,Wabef,Wmnij,WmBeJ,WmBEj,oovv,ovov,ovvv,vvov,ooov,tia,tiatia,tijab,Dijab)
    _tia = permutedims(tia,[2,1])
    @tensor begin
        _tijab[i,j,a,b] := (oovv[i,j,a,b] 
                           + tijab[i,j,a,e]*(Fae[b,e] - 0.5*tia[m,b]*Fme[m,e])
                           + tijab[i,j,e,b]*(Fae[a,e] - 0.5*tia[m,a]*Fme[m,e])
                           - tijab[i,m,a,b]*(Fmi[m,j] + 0.5*tia[j,e]*Fme[m,e])
                           - tijab[m,j,a,b]*(Fmi[m,i] + 0.5*tia[i,e]*Fme[m,e])
                           + (tijab[m,n,a,b] + tiatia[m,n,a,b])*Wmnij[m,n,i,j]
                           + (tijab[i,j,e,f] + tiatia[i,j,e,f])*Wabef[a,b,e,f]
                           + ((tijab[i,m,a,e] - tijab[m,i,a,e])*WmBeJ[m,b,e,j]
                              - tiatia[i,m,e,a]*oovv[m,j,e,b])
                           + tijab[i,m,a,e]*(WmBeJ[m,b,e,j] + WmBEj[m,b,e,j])
                           + (tijab[m,i,b,e]*WmBEj[m,a,e,j] 
                              - tiatia[i,m,e,b]*ovov[m,a,j,e])
                           + (tijab[m,j,a,e]*WmBEj[m,b,e,i]
                              - tiatia[j,m,e,a]*ovov[m,b,i,e])
                           + ((tijab[j,m,b,e] - tijab[m,j,b,e])*WmBeJ[m,a,e,i]
                              - tiatia[j,m,e,b]*oovv[m,i,e,a])
                           + tijab[j,m,b,e]*(WmBeJ[m,a,e,i] + WmBEj[m,a,e,i])
                           + _tia[e,i]*ovvv[j,a,b,e]
                           + _tia[e,j]*ovvv[i,b,a,e]
                           - tia[m,a]*ooov[m,j,i,b]
                           - tia[m,b]*ooov[j,i,m,a]
                          )
    end
    _tijab .= _tijab ./ Dijab
    return _tijab
end
end #module
