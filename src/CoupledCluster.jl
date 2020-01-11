module CoupledCluster
"""
Basic module for running CC computations in Julia.

short term goal is RCCD and UCCD
medium term goal is RCCSD and UCCSD
long term goal is RCCSD(T) and UCCSD(T)

Implemented --> RCCD
			--> 
Optimized --> RCCD


usage --> methods should be defined like do_<r/u><method> and take in
	  --> a Wavefunction.jl Wfn object as their sole _required_ input.
	  --> optional inputs such as maxit, convergence, etc can be defined
	  --> via multiple dispatch
"""

using JuES.Wavefunction
using JuES.Transformation
using Base.Threads
using LinearAlgebra
using Dates
export do_rccd


function do_rccd(refWfn::Wfn)
    #implicit maxit = 40
    return do_rccd(refWfn, 40)
end
@fastmath @inbounds function do_rccd(refWfn::Wfn, maxit, doprint = false)
    #goes through appropriate steps to do RCCD
    set_zero_subnormals(true)
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    ovov = refWfn.ijab
    vvvv =
        tei_transform(refWfn.uvsr, refWfn.Cav, refWfn.Cav, refWfn.Cav, refWfn.Cav, "test")
    ovvo =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cav, refWfn.Cav, refWfn.Cao, "test")
    oooo =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cao, refWfn.Cao, "test")
    ooov =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cao, refWfn.Cav, "test")
    oovv =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cav, refWfn.Cav, "test")
    dtt = eltype(ovov)
    epsa = refWfn.epsa
    T2 = zeros(dtt, nocc, nocc, nvir, nvir)
    Dijab = form_Dijab(T2, epsa)
    T2_init!(T2, ovov, Dijab)
    println("@RMP2 ", ccenergy(T2, ovov))
    Fae = form_Fae(T2, ovov)
    Fmi = form_Fmi(T2, ovov)
    Wmnij = form_Wmnij(oooo, ovov, T2)
    Wmnij = zeros(size(Wmnij))
    Wabef = form_Wabef(vvvv, ovov, T2)
    Wabef = zeros(size(Wabef))
    WmBeJ = form_WmBeJ(ovvo, ovov, T2)
    WmBeJ = zeros(size(WmBeJ))
    WmBEj = form_WmBEj(ovov, oovv, T2)
    WmBEj = zeros(size(WmBEj))
    dt = @elapsed for i in UnitRange(1, maxit) #TODO: implement RMS check
        T2 = cciter(
            T2,
            ovov,
            vvvv,
            oooo,
            oovv,
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
            println("@CCD ", ccenergy(T2, ovov))
        end
        t1 = Dates.Time(Dates.now())
    end
    if doprint
        println("CCD energy computed in $dt s")
    end
    return ccenergy(T2, ovov)
end

function ccenergy(tiJaB, iajb)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for i in rocc
        for j in rocc
            for a in rvir
                for b in rvir
                    cache = iajb[i, a, j, b]
                    ecc += cache * 2 * tiJaB[i, j, a, b]
                    ecc -= cache * tiJaB[j, i, a, b]
                end
            end
        end
    end
    return ecc
end

function cciter(
    tiJaB_i,
    ovov,
    vvvv,
    oooo,
    oovv,
    ovvo,
    Dijab,
    Fae,
    Fmi,
    Wabef,
    Wmnij,
    WmBeJ,
    WmBEj,
)
    form_Fae!(Fae, tiJaB_i, ovov)
    #Fae = zeros(size(Fae))
    form_Fmi!(Fmi, tiJaB_i, ovov)
    #Fmi = zeros(size(Fmi))
    form_Wmnij!(Wmnij, oooo, ovov, tiJaB_i)
    #Wmnij = zeros(size(Wmnij))
    form_Wabef!(Wabef, vvvv, ovov, tiJaB_i)
    #Wabef = zeros(size(Wabef))
    form_WmBeJ!(WmBeJ, ovvo, ovov, tiJaB_i)
    #WmBeJ = zeros(size(WmBeJ))
    WmBEj = form_WmBEj!(WmBEj, ovov, oovv, tiJaB_i)
    #WmBEj = zeros(size(WmBEj))
    tiJaB_d = form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, ovov, Dijab)
    return tiJaB_d
end

function T2_init!(tiJaB, iajb, Dijab)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    tiJaB .= permutedims(iajb, [1, 3, 2, 4]) ./ Dijab
end

function form_Fae(tiJaB, menf)
    dt = eltype(menf)
    nvir = size(tiJaB, 4)
    Fae = zeros(dt, nvir, nvir)
    form_Fae!(Fae, tiJaB, menf)
    return Fae
end
function form_Fae!(Fae, tiJaB, menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    Fae .= 0.0
    cache1 = zeros(eltype(menf), nocc, nocc)
    cache2 = zeros(eltype(menf), nocc, nocc)
    for f in rvir
        for a in rvir
            for e in rvir
                for n in rocc
                    @simd for m in rocc
                        Fae[a, e] -=
                            menf[m, e, n, f] * (2 * tiJaB[m, n, a, f] - tiJaB[n, m, a, f])
                    end
                end
            end
        end
    end
end

function form_Fmi(tiJaB, menf)
    dt = eltype(menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    Fmi = zeros(dt, nocc, nocc)
    form_Fmi!(Fmi, tiJaB, menf)
    return Fmi
end
function form_Fmi!(Fmi, tiJaB, menf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    Fmi .= 0.0
    for f in rvir
        for e in rvir
            for n in rocc
                for i in rocc
                    @simd for m in rocc
                        Fmi[m, i] += menf[m, e, n, f] 
                            * (2 * tiJaB[i, n, e, f] - tiJaB[i, n, f, e])
                    end
                end
            end
        end
    end
    #return Fmi
end

function form_Dijab(tiJaB, F)
    dt = eltype(tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    Dijab = zeros(dt, nocc, nocc, nvir, nvir)
    for i in rocc
        for j in rocc
            for a in rvir
                for b in rvir
                    aa = a + nocc
                    bb = b + nocc
                    Dijab[i, j, a, b] = F[i] + F[j] - F[aa] - F[bb]
                end
            end
        end
    end
    return Dijab
end

function form_T2(tiJaB_i, Fae, Fmi, WmBeJ, WmBEj, Wabef, Wmnij, iajb, Dijab)
    dtt = eltype(tiJaB_i)
    nocc = size(Wmnij, 1)
    nvir = size(tiJaB_i, 4)
    tiJaB_d = zeros(dtt, nocc, nocc, nvir, nvir)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for b in rvir
        for a in rvir
            for j in rocc
                for i in rocc
                    temp = iajb[i, a, j, b]
                    for e in rvir
                        temp += tiJaB_i[i, j, a, e] * Fae[b, e]
                        temp += tiJaB_i[j, i, b, e] * Fae[a, e]
                        for f in rvir
                            temp += tiJaB_i[i, j, e, f] * Wabef[a, b, e, f]
                        end
                    end
                    for m in rocc
                        temp -= tiJaB_i[i, m, a, b] * Fmi[m, j]
                        temp -= tiJaB_i[m, j, a, b] * Fmi[m, i]
                        for n in rocc
                            temp += tiJaB_i[m, n, a, b] * Wmnij[m, n, i, j]
                        end
                        for e in rvir
                            temp += 2 * tiJaB_i[i, m, a, e] * WmBeJ[m, b, e, j]
                            temp -= tiJaB_i[m, i, a, e] * WmBeJ[m, b, e, j]
                            temp += tiJaB_i[i, m, a, e] * WmBEj[m, b, e, j]
                            temp += tiJaB_i[m, i, b, e] * WmBEj[m, a, e, j]
                            temp += tiJaB_i[m, j, a, e] * WmBEj[m, b, e, i]
                            temp += 2 * tiJaB_i[j, m, b, e] * WmBeJ[m, a, e, i]
                            temp -= tiJaB_i[m, j, b, e] * WmBeJ[m, a, e, i]
                            temp += tiJaB_i[j, m, b, e] * WmBEj[m, a, e, i]
                        end
                    end
                    tiJaB_d[i, j, a, b] = temp
                end
            end
        end
    end
    tiJaB_d .= tiJaB_d ./ Dijab
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
function form_Wmnij!(Wmnij, minj, menf, tiJaB)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for j in rocc
        for n in rocc
            for i in rocc
                for m in rocc
                    Wmnij[m, n, i, j] = minj[m, i, n, j]
                    for f in rvir
                        for e in rvir
                            Wmnij[m, n, i, j] += tiJaB[i, j, e, f] * menf[m, e, n, f] / 2.0
                        end
                    end
                end
            end
        end
    end
end

function form_Wabef(aebf, mnef, tiJaB)
    dt = eltype(aebf)
    nvir = size(tiJaB, 4)
    Wabef = zeros(dt, nvir, nvir, nvir, nvir)
    form_Wabef!(Wabef, aebf, mnef, tiJaB)
    return Wabef
end
function form_Wabef!(Wabef, aebf, menf, tiJaB)
    dtt = eltype(Wabef)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    _iJaB = zeros(dtt, nocc, nocc)
    Wabef .= 0
    for f in rvir
        for e in rvir
            for b in rvir
                for a in rvir
                    Wabef[a, b, e, f] = aebf[a, e, b, f]
                    for n in rocc
                        for m in rocc
                            Wabef[a, b, e, f] += tiJaB[m, n, a, b] * menf[m, e, n, f] / 2.0
                        end
                    end
                end
            end
        end
    end
end

function form_WmBeJ(mebj, iajb, tiJaB)
    dtt = eltype(iajb)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBeJ = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBeJ!(WmBeJ, mebj, iajb, tiJaB)
    return WmBeJ
end
function form_WmBeJ!(WmBeJ, mebj, iajb, tiJaB)
    dtt = eltype(WmBeJ)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    #WmBeJ .= 0 
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    _iJaB = zeros(dtt, nocc, nocc)
    for e in rvir
        for b in rvir
            for j in rocc
                for m in rocc
                    WmBeJ[m, b, e, j] = mebj[m, e, b, j]
                    for f in rvir
                        for n in rocc
                            WmBeJ[m, b, e, j] +=
                                0.5 * iajb[m, e, n, f] *
                                (2 * tiJaB[n, j, f, b] - tiJaB[j, n, f, b]) / 2.0
                            WmBeJ[m, b, e, j] -= iajb[n, e, m, f] * tiJaB[n, j, f, b] / 2.0
                        end
                    end
                end
            end
        end
    end
end


function form_WmBEj(nemf, mjbe, tiJaB)
    dtt = eltype(nemf)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    WmBEj = zeros(dtt, nocc, nvir, nvir, nocc)
    form_WmBEj!(WmBEj, nemf, mjbe, tiJaB)
    return WmBEj
end
function form_WmBEj!(WmBEj, nemf, mjbe, tiJaB)
    dtt = eltype(WmBEj)
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for e in rvir
        for m in rocc
            for b in rvir
                for j in rocc
                    WmBEj[m, b, e, j] = -mjbe[m, j, b, e]
                    for f in rvir
                        for n in rocc
                            WmBEj[m, b, e, j] += tiJaB[j, n, f, b] * nemf[n, e, m, f] / 2.0
                        end
                    end
                end
            end
        end
    end
    return WmBEj
end
end #module CC
