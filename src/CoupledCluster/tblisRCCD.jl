"""
    JuES.CoupledCluster.RCCD

performs coupled cluster doubles (CCD) computations.

## Functions

    JuES.CoupledCluster.RCCD.do_rccd
"""
module tRCCD
using JuES.Output
using JuES.Wavefunction
using TensorOperations
using TBLIS
using JuES.IntegralTransformation
using JuES
include("Denominators.jl")
export do_trccd
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
function do_trccd(refWfn::Wfn; kwargs...)
    maxit = 40
    doprint = false
    return_T2 = false
    TBLIS.init()
    JuES.CoupledCluster.print_header()
    t = @elapsed begin
        nocc = refWfn.nalpha
        nvir = refWfn.nvira
        epsa = refWfn.epsa
        T = eltype(refWfn.ao_eri)
        oovv,ovov,ovvo,oooo,vvvv = make_rccd_integrals(refWfn,refWfn.Cao,refWfn.Cav) 
        T2 = zeros(T, nocc, nocc, nvir, nvir)
        Dijab = form_Dijab(T2, epsa)
        #T2_init!(T2, ovov, Dijab)
        #@output "RMP2 {:>20.17f}\n" ccenergy(T2,oovv)
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
            @output "{:1} {:>20.17f}\n" i ccenergy(T2,oovv)
        end
        if doprint
            println("CCD iterations computed in $dt s")
        end
    end
    @output "CCD completed in {:>5.2f}s\n" t
    if return_T2
        return ccenergy(T2, oovv), T2
    else
        return ccenergy(T2, oovv)
    end
end

function make_rccd_integrals(wfn,Cao,Cav)
    oovv = get_eri(wfn,"OOVV")
    vvvv = get_eri(wfn,"VVVV")
    ovvo = get_eri(wfn,"OVVO")
    ovov = get_eri(wfn,"OVOV")
    oooo = get_eri(wfn,"OOOO") 
    #oovv = tei_transform(gao, Cao, Cav, Cao, Cav, "oovv")
    #vvvv = tei_transform(gao, Cav, Cav, Cav, Cav, "vvvv")
    #ovvo = tei_transform(gao, Cao, Cav, Cav, Cao, "ovvo")
    #ovov = tei_transform(gao, Cao, Cao, Cav, Cav, "oovv")
    #oooo = tei_transform(gao, Cao, Cao, Cao, Cao, "oooo")
    #oovv = permutedims(oovv,[1,3,2,4])
    #ovov = permutedims(ovov,[1,3,2,4])
    #ovvo = permutedims(ovvo,[1,3,2,4])
    #oooo = permutedims(oooo,[1,3,2,4])
    #vvvv = permutedims(vvvv,[1,3,2,4])
    return oovv,ovov,ovvo,oooo,vvvv
end

function ccenergy(tiJaB, ijab)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for i in rocc
        for j in rocc
            #cache = ijab[i,j,:,:]
            for a in rvir
                for b in rvir
                    ecc += ijab[i,j,a,b] * 2 * tiJaB[i, j, a, b]
                    ecc -= ijab[i,j,a,b] * tiJaB[j, i, a, b]
                end
            end
        end
    end
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
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    for i in rocc
        for j in rocc
            cache = iajb[i,:,j,:]
            for a in rvir
                for b in rvir
                    tiJaB[i,j,a,b] = cache[a,b] / Dijab[i,j,a,b]
                end
            end
        end
    end
end

function form_Fae(tiJaB, menf)
    dt = eltype(menf)
    nvir = size(tiJaB, 4)
    Fae = zeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, menf)
    return Fae
end
function form_Fae!(Fae, tiJaB, mnef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = UnitRange(1, nocc)
    rvir = UnitRange(1, nvir)
    Fae .= 0.0
    temp = zeros(eltype(tiJaB),size(tiJaB))
    _temp = TTensor{eltype(tiJaB)}(temp)
    _tijab1 = TTensor{eltype(tiJaB)}(tiJaB,2.0)
    _tijab2 = TTensor{eltype(tiJaB)}(tiJaB,-1.0)
    _mnef = TTensor{eltype(tiJaB)}(mnef,-1.0)
    _Fae = TTensor{eltype(tiJaB)}(Fae)
    add!(_tijab1,_temp,"mnaf","mnaf")
    add!(_tijab2,_temp,"nmaf","mnaf")
    mul!(_Fae,_mnef,_temp,"mnef","mnaf","ae")
    #@tensoropt begin
    #    Fae[a,e] = -1*(mnef[m,n,e,f]*(2*tiJaB[m,n,a,f] - tiJaB[n,m,a,f]))
    #end
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
    Fmi .= 0.0
    @tensoropt begin
        Fmi[m,i] = mnef[m,n,e,f]*(2*tiJaB[i,n,e,f] - tiJaB[i,n,f,e])
    end
    return Fmi
end

function form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, ijab, Dijab)
    dtt = eltype(tiJaB_i)
    nocc = size(Wmnij, 1)
    nvir = size(tiJaB_i, 4)
    rocc = 1:nocc
    rvir = 1:nvir 
    #ijab = permutedims(iajb,[1,3,2,4])
    tiJaB_d = deepcopy(ijab)
    _tiJaB_d = TTensor{eltype(tiJaB_i)}(tiJaB_d,1.0)
    _tiJaB_i = TTensor{eltype(tiJaB_i)}(tiJaB_i,1.0)
    _Fae = TTensor{eltype(tiJaB_i)}(Fae)
    _Fmi = TTensor{eltype(tiJaB_i)}(Fmi,-1.0)
    _Wmnij = TTensor{eltype(tiJaB_i)}(Wmnij)
    _Wabef = TTensor{eltype(tiJaB_i)}(Wabef)
    _WmBeJ = TTensor{eltype(tiJaB_i)}(WmBeJ)
    _2WmBeJ = TTensor{eltype(tiJaB_i)}(WmBeJ,2.0)
    _n1WmBeJ = TTensor{eltype(tiJaB_i)}(WmBeJ,-1.0)
    _WmBEj = TTensor{eltype(tiJaB_i)}(WmBEj)

    mul!(_tiJaB_d,_tiJaB_i,_Fae,"ijae","be","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_Fae,"jibe","ae","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_Fmi,"imab","mj","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_Fmi,"mjab","mi","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_Wmnij,"mnab","mnij","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_Wabef,"ijef","abef","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_2WmBeJ,"imae","mbej","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_n1WmBeJ,"miae","mbej","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_WmBEj,"imae","mbej","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_WmBEj,"mibe","maej","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_WmBEj,"mjae","mbei","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_2WmBeJ,"jmbe","maei","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_n1WmBeJ,"mjbe","maei","ijab")
    mul!(_tiJaB_d,_tiJaB_i,_WmBEj,"jmbe","maei","ijab")
    tiJaB_d = tiJaB_d ./ Dijab
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
    _Wmnij = TTensor{eltype(tiJaB)}(Wmnij)
    _p5mnef = TTensor{eltype(tiJaB)}(mnef,0.5)
    _tiJaB = TTensor{eltype(tiJaB)}(tiJaB)
    mul!(_Wmnij,_tiJaB,_p5mnef,"ijef","mnef","mnij")
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
    _Wabef = TTensor{eltype(tiJaB)}(Wabef)
    _tiJaB = TTensor{eltype(tiJaB)}(tiJaB)
    _p5mnef = TTensor{eltype(tiJaB)}(mnef,0.5)
    mul!(_Wabef,_tiJaB,_p5mnef,"mnab","mnef","abef")
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
    _WmBeJ = TTensor{eltype(tiJaB)}(WmBeJ)
    _temp = TTensor{eltype(tiJaB)}(zeros(eltype(tiJaB),size(tiJaB)))
    _2tiJaB = TTensor{eltype(tiJaB)}(tiJaB,2.0)
    _tiJaB = TTensor{eltype(tiJaB)}(tiJaB,1.0)
    _ntiJaB = TTensor{eltype(tiJaB)}(tiJaB,-1.0)
    _p5mnef = TTensor{eltype(tiJaB)}(mnef,0.5)
    _np5mnef = TTensor{eltype(tiJaB)}(mnef,-0.5)
    add!(_2tiJaB,_temp,"njfb","njfb")
    add!(_ntiJaB,_temp,"jnfb","njfb")
    mul!(_WmBeJ,_p5mnef,_temp,"mnef","njfb","mbej")
    mul!(_WmBeJ,_np5mnef,_tiJaB,"nmef","njfb","mbej")
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
    WmBEj .= 0
    _mbje = TTensor{eltype(tiJaB)}(mbje,-1.0)
    _WmBEj = TTensor{eltype(tiJaB)}(WmBEj)
    add!(_mbje,_WmBEj,"mbje","mbej")
    _p5tiJaB = TTensor{eltype(tiJaB)}(tiJaB,0.5)
    _nmef = TTensor{eltype(tiJaB)}(nmef)
    mul!(_WmBEj,_p5tiJaB,_nmef,"jnfb","nmef","mbej")
    return WmBEj
end
end #module
