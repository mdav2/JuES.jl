"""
    JuES.CoupledCluster.mRCCD

performs coupled cluster doubles (CCD) computations using explicit matrix multiplications.

## Functions

    JuES.CoupledCluster.RCCD.do_rccd
"""
module mRCCD
using Base.Threads
using JuES.Wavefunction
using CuArrays
using JuES
using TensorOperations
using LinearAlgebra
#BLAS.set_num_threads(12)
#using TensorCast
using JuES.Transformation
include("Denominators.jl")
#include("GSGEMM.jl")
export do_rccd
"""
    do_rccd

Applies RCCD equations to the input Wfn object. Disk based versus in-core 
algorithm is selected based on the type of atomic orbitals in Wfn.uvsr.
---
## paramters
refWfn::Wfn         -> Wfn object to which the RCCD equations will be applied.

maxit::Int          -> maximum number of coupled cluster iterations.

doprint::Bool=false -> whether or not to print energy and timing information to 
    stdout.
## output
ccenergy::Float -> final RCCD energy. 
"""
function do_rccd(refWfn::Wfn; maxit=40, doprint=false, return_T2=false)
    JuES.CoupledCluster.print_header()
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    epsa = refWfn.epsa
    T = eltype(refWfn.ao_eri)
    oovv,ovov,ovvo,oooo,vvvv = make_rccd_integrals(refWfn.ao_eri,refWfn.Cao,refWfn.Cav) 
    T2 = zeros(T, nocc, nocc, nvir, nvir)
    Dijab = form_Dijab(T2, epsa)
    T2_init!(T2, ovov, Dijab)
    if doprint
        println("@RMP2 ", ccenergy(T2, oovv))
    end
    oovv = Array{T}(oovv)
    ovov = Array{T}(ovov)
    ovvo = Array{T}(ovvo)
    oooo = Array{T}(oooo)
    vvvv = Array{T}(vvvv)
    Fae = form_Fae(T2, oovv)
    Fmi = form_Fmi(T2, oovv)
    Wmnij = form_Wmnij(oooo, oovv, T2)
    Wabef = form_Wabef(vvvv, oovv, T2)
    WmBeJ = form_WmBeJ(ovvo, oovv, T2)
    WmBEj = form_WmBEj(oovv, ovov, T2)
    dt = @elapsed for i in 0:maxit-1 #TODO: implement RMS check
        T2 = cciter(
            T2,
            oovv,
            vvvv,
            oooo,
            ovov,
            ovvo,
            Dijab,
            Fae,
            Fmi,
            Wabef,
            Wmnij,
            WmBeJ,
            WmBEj,
        )
        if doprint
            println("$i @CCD ", ccenergy(T2, oovv))
        end
    end
    if doprint
        println("CCD iterations computed in $dt s")
    end
    if return_T2
        return ccenergy(T2, oovv), T2
    else
        return ccenergy(T2, oovv)
    end
end

function make_rccd_integrals(gao::Array,Cao,Cav)
    oovv = tei_transform(gao, Cao, Cav, Cao, Cav, "oovv")
    vvvv = tei_transform(gao, Cav, Cav, Cav, Cav, "vvvv")
    ovvo = tei_transform(gao, Cao, Cav, Cav, Cao, "ovvo")
    ovov = tei_transform(gao, Cao, Cao, Cav, Cav, "oovv")
    oooo = tei_transform(gao, Cao, Cao, Cao, Cao, "oooo")
    oovv = permutedims(oovv,[1,3,2,4])
    ovov = permutedims(ovov,[1,3,2,4])
    ovvo = permutedims(ovvo,[1,3,2,4])
    oooo = permutedims(oooo,[1,3,2,4])
    vvvv = permutedims(vvvv,[1,3,2,4])
    return oovv,ovov,ovvo,oooo,vvvv
end

function ccenergy(tiJaB, ijab)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    @tensor ecc = ijab[i,j,a,b]*2*tiJaB[i,j,a,b] - ijab[i,j,a,b]*tiJaB[j,i,a,b]
    return ecc
end

@inbounds @fastmath function cciter(
    tiJaB_i,
    oovv,
    vvvv,
    oooo,
    ovov,
    ovvo,
    Dijab,
    Fae,
    Fmi,
    Wabef,
    Wmnij,
    WmBeJ,
    WmBEj,
)
    form_Wabef!(Wabef, vvvv, oovv, tiJaB_i)
    form_WmBeJ!(WmBeJ, ovvo, oovv, tiJaB_i)
    form_WmBEj!(WmBEj, oovv, ovov, tiJaB_i)
    form_Wmnij!(Wmnij, oooo, oovv, tiJaB_i)
    form_Fae!(Fae, tiJaB_i, oovv)
    form_Fmi!(Fmi, tiJaB_i, oovv)

    tiJaB_d = form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, oovv, Dijab)
    return tiJaB_d
end

function T2_init!(tiJaB, iajb, Dijab)
    #nocc = size(tiJaB, 1)
    #nvir = size(tiJaB, 4)
    #rocc = UnitRange(1, nocc)
    #rvir = UnitRange(1, nvir)
    tiJaB .= permutedims(iajb,(1,3,2,4)) ./ Dijab
    #for i in rocc
    #    for j in rocc
    #        cache = iajb[i,:,j,:]
    #        for a in rvir
    #            for b in rvir
    #                tiJaB[i,j,a,b] = cache[a,b] / Dijab[i,j,a,b]
    #            end
    #        end
    #    end
    #end
end

function form_Fae(tiJaB, menf)
    dt = eltype(menf)
    nvir = size(tiJaB, 4)
    Fae = zeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, menf)
    return Fae
end
function form_Fae!(Fae, tiJaB, mnef)
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _tiJaB = 2*reshape(permutedims(tiJaB,(3,1,2,4)),v,o^2*v)
    _tiJaB -= reshape(permutedims(tiJaB,(3,2,1,4)),v,o^2*v)
    _mnef = reshape(permutedims(mnef,(1,2,4,3)),o^2*v,v)
    BLAS.gemm!('N','N',-1.0f0,_tiJaB,_mnef,0.0f0,Fae)
end

function form_Fmi(tiJaB, menf)
    dt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Fmi = zeros(dt, nocc, nocc)
    form_Fmi!(Fmi, tiJaB, menf)
    return Fmi
end
function form_Fmi!(Fmi, tiJaB, mnef)
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    #@cast _tiJaB[(n,e,f),(i)] := 2*tiJaB[i,n,e,f] - tiJaB[i,n,f,e]
    _tiJaB = 2*reshape(permutedims(tiJaB,(2,3,4,1)),o*v^2,o)
    _tiJaB -= reshape(permutedims(tiJaB,(2,4,3,1)),o*v^2,o)
    #@cast _mnef[(m),(n,e,f)] := mnef[m,n,e,f]
    _mnef = reshape(mnef,o,o*v^2)
    #mul!(Fmi,_mnef,_tiJaB)
    BLAS.gemm!('N','N',1.0f0,_mnef,_tiJaB,0.0f0,Fmi)
    return Fmi
end

@fastmath @inbounds function form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, ijab, Dijab)
    chonk = 1
    T = eltype(tiJaB_i)
    tiJaB_d = zeros(T,size(ijab[:,:,:,:]))
    tiJaB_d .= ijab
    o = size(ijab,1)
    v = size(ijab,3)
    _tiJaB_i = reshape(tiJaB_i,o^2*v,v)
    _tiJaB_d = reshape(tiJaB_d,o^2*v,v)
    BLAS.gemm!('N','T',1.0f0,_tiJaB_i,Fae,1.0f0,_tiJaB_d)
    tiJaB_d .= reshape(_tiJaB_d,o,o,v,v)

    @tensor _tiJaB_d[j,i,b,a] := tiJaB_d[i,j,a,b]
    #_tiJaB_d = reshape(permutedims(tiJaB_d,(2,1,4,3)),o^2*v,v)
    _tiJaB_d = reshape(_tiJaB_d,o^2*v,v)
    _tiJaB_i = reshape(tiJaB_i,o^2*v,v)
    BLAS.gemm!('N','T',1.0f0,_tiJaB_i,Fae,1.0f0,_tiJaB_d)
    tiJaB_d .= permutedims(reshape(_tiJaB_d,o,o,v,v),(2,1,4,3))

    @tensor _tiJaB_d[1,3,4,2] := tiJaB_d[1,2,3,4]
    _tiJaB_d = reshape(_tiJaB_d,o*v^2,o)
    _tiJaB_i = reshape(permutedims(tiJaB_i,(2,1,3,4)),o,o*v^2)
    BLAS.gemm!('T','N',-1.0f0,_tiJaB_i,Fmi,1.0f0,_tiJaB_d)
    tiJaB_d .= permutedims(reshape(transpose(_tiJaB_d),o,o,v,v),(2,1,3,4))

    #_tiJaB_d = reshape(permutedims(tiJaB_d,(2,3,4,1)),o*v^2,o)
    @tensor _tiJaB_d2[2,3,4,1] := tiJaB_d[1,2,3,4]
    _tiJaB_d2 = reshape(_tiJaB_d2,o*v^2,o)
    _tiJaB_i2 = reshape(tiJaB_i,o,o*v^2)
    BLAS.gemm!('T','N',-1.0f0,_tiJaB_i2,Fmi,1.0f0,_tiJaB_d2)
    tiJaB_d .= reshape(transpose(_tiJaB_d2),o,o,v,v)

    _tiJaB_d = reshape(tiJaB_d,o^2,v^2)
    _tiJaB_i = reshape(tiJaB_i,o^2,v^2)
    _Wmnij = reshape(Wmnij,o^2,o^2)
    BLAS.gemm!('T','N',1.0f0,_Wmnij,_tiJaB_i,1.0f0,_tiJaB_d)

    _Wabef = reshape(Wabef,v^2,v^2)
    BLAS.gemm!('N','T',1.0f0,_tiJaB_i,_Wabef,1.0f0,_tiJaB_d)
    #GSGEMM!(_tiJaB_d,1.0f0,_tiJaB_i,_Wabef,1.0f0,size(_tiJaB_i,2),size(_tiJaB_i,2))
    tiJaB_d = reshape(_tiJaB_d,o,o,v,v)

    scr1 = Array{T}(undef,o*v,o*v)
    scr2 = Array{T}(undef,o*v,o*v)
    scr3 = Array{T}(undef,o*v,o*v)

    @tensor _tiJaB_d[1,3,4,2] := tiJaB_d[1,2,3,4]
    _tiJaB_d = reshape(_tiJaB_d,o*v,v*o)
    @tensor _scr2[1,3,2,4] := tiJaB_i[1,2,3,4]
    scr2 = reshape(_scr2,o*v,o*v)
    @tensor _scr3[1,3,2,4] := WmBeJ[1,2,3,4]
    scr3 = 2*reshape(_scr3,o*v,v*o)

    @tensor int[1,3,2,4] := WmBEj[1,2,3,4]
    int = reshape(int,o*v,v*o)
    scr3 += int
    BLAS.gemm!('N','N',1.0f0,scr2,scr3,1.0f0,_tiJaB_d)

    scr3 -= int
    scr3 *= 0.5f0
    scr2 = reshape(permutedims(tiJaB_i,(1,4,2,3)),o*v,o*v)
    BLAS.gemm!('T','N',-1.0f0,scr2,scr3,1.0f0,_tiJaB_d)

    tiJaB_d = permutedims(reshape(_tiJaB_d,o,v,v,o),(1,4,2,3))

    _tiJaB_d = reshape(permutedims(tiJaB_d,(1,4,3,2)),o*v,v*o)
    scr2 = reshape(permutedims(tiJaB_i,(1,4,2,3)),o*v,o*v)
    scr4 = reshape(permutedims(WmBEj,(1,3,2,4)),o*v,o*v)
    tmp = Array{T}(undef,size(_tiJaB_d))
    BLAS.gemm!('T','N',1.0f0,scr2,scr4,0.0f0,tmp)
    _tiJaB_d += tmp
    tiJaB_d = permutedims(reshape(_tiJaB_d,o,v,v,o),(1,4,3,2))
    

    scr1 = reshape(permutedims(tiJaB_d,(2,3,4,1)),o*v,v*o)
    scr1 += tmp
    tiJaB_d = permutedims(reshape(scr1,o,v,v,o),(4,1,2,3))

    scr1 = reshape(permutedims(tiJaB_d,(2,4,3,1)),o*v,v*o)
    scr2 = reshape(permutedims(tiJaB_i,(1,3,2,4)),o*v,o*v)
    BLAS.gemm!('N','N',1.0f0,scr2,scr4,1.0f0,scr1)

    scr2 *= 2
    scr2 -= reshape(permutedims(tiJaB_i,(2,3,1,4)),o*v,o*v)
    BLAS.gemm!('N','N',1.0f0,scr2,scr3,1.0f0,scr1)

    tiJaB_d = permutedims(reshape(scr1,o,v,v,o),(4,1,3,2)) ./ Dijab

    return tiJaB_d
end
function form_Wmnij(minj, menf, tiJaB)
    dtt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Wmnij = zeros(dtt, nocc, nocc, nocc, nocc)
    form_Wmnij!(Wmnij, minj, menf, tiJaB)
    return Wmnij
end
function form_Wmnij!(Wmnij, mnij, mnef, tiJaB)
    Wmnij .= mnij
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _Wmnij = reshape(Wmnij,o^2,o^2)
    _tiJaB = reshape(tiJaB,o^2,v^2)
    _mnef = reshape(mnef,o^2,v^2)
    BLAS.gemm!('N','T',0.5f0,_mnef,_tiJaB,1.0f0,_Wmnij)
    Wmnij .= reshape(_Wmnij,o,o,o,o)
end

function form_Wabef(aebf, mnef, tiJaB)
    dt = eltype(aebf)
    nvir = size(tiJaB, 4)
    Wabef = zeros(dt, nvir, nvir, nvir, nvir)
    form_Wabef!(Wabef, aebf, mnef, tiJaB)
    return Wabef
end
function form_Wabef!(Wabef, abef, mnef, tiJaB)
    Wabef .= abef
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _Wabef = reshape(Wabef,v^2,v^2)
    _tiJaB = reshape(tiJaB,o^2,v^2)
    _mnef = reshape(mnef,o^2,v^2)
    BLAS.gemm!('T','N',0.5f0,_tiJaB,_mnef,1.0f0,_Wabef)
    Wabef = reshape(_Wabef,v,v,v,v)
    #@tensoropt begin
    #    Wabef[a,b,e,f] = tiJaB[m,n,a,b]*mnef[m,n,e,f]/2 + abef[a,b,e,f]
    #end
end

function form_WmBeJ(mebj, iajb, tiJaB)
    dtt = eltype(iajb)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBeJ = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBeJ!(WmBeJ, mebj, iajb, tiJaB)
    return WmBeJ
end
function form_WmBeJ!(WmBeJ, mbej, mnef, tiJaB)
    WmBeJ .= mbej
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _WmBeJ = reshape(permutedims(WmBeJ,(1,3,4,2)),o*v,o*v)
    _tiJaB1 = 2*reshape(permutedims(tiJaB,(1,3,2,4)),o*v,o*v)
    _tiJaB1 -= reshape(permutedims(tiJaB,(2,3,1,4)),o*v,o*v)
    _mnef1 = reshape(permutedims(mnef,(1,3,2,4)),o*v,o*v)

    _mnef2 = reshape(permutedims(mnef,(1,4,2,3)),o*v,o*v)
    _tiJaB2 = reshape(permutedims(tiJaB,(1,3,2,4)),o*v,o*v)
    BLAS.gemm!('N','N',0.5f0,_mnef1,_tiJaB1,1.0f0,_WmBeJ)
    BLAS.gemm!('T','N',-0.5f0,_mnef2,_tiJaB2,1.0f0,_WmBeJ)
    WmBeJ .= permutedims(reshape(_WmBeJ,o,v,o,v),(1,4,2,3))
end


function form_WmBEj(nemf, mjbe, tiJaB)
    dtt = eltype(nemf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBEj = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBEj!(WmBEj, nemf, mjbe, tiJaB)
    return WmBEj
end
function form_WmBEj!(WmBEj, nmef, mbje, tiJaB)
    WmBEj .= -permutedims(mbje,(1,2,4,3))
    o = size(tiJaB,1)
    v = size(tiJaB,3)
    _tiJaB = reshape(permutedims(tiJaB,(2,3,1,4)),o*v,o*v)
    _nmef = reshape(permutedims(nmef,(2,3,1,4)),o*v,o*v)
    _WmBEj = (reshape(permutedims(WmBEj,(1,3,4,2)),o*v,o*v))
    BLAS.gemm!('N','N',0.5f0,_nmef,_tiJaB,1.0f0,_WmBEj)
    WmBEj .= permutedims(reshape(_WmBEj,o,v,o,v),(1,4,2,3))
    return WmBEj
end
end #module
