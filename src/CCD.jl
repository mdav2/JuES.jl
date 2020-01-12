function do_rccd(refWfn::Wfn)
    #implicit maxit = 40
    return do_rccd(refWfn, 40)
end
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
function do_rccd(refWfn::Wfn, maxit; doprint::Bool=false, return_T2::Bool=false)
    #goes through appropriate steps to do RCCD
    set_zero_subnormals(true)
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    ovov = refWfn.ijab
    vvvv =
        tei_transform(refWfn.uvsr, refWfn.Cav, refWfn.Cav, refWfn.Cav, refWfn.Cav, "vvvv")
    ovvo =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cav, refWfn.Cav, refWfn.Cao, "ovvo")
    oooo =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cao, refWfn.Cao, "oooo")
    ooov =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cao, refWfn.Cav, "ooov")
    oovv =
        tei_transform(refWfn.uvsr, refWfn.Cao, refWfn.Cao, refWfn.Cav, refWfn.Cav, "oovv")
    dtt = Float64#eltype(ovov)
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
    if return_T2
        return ccenergy(T2, ovov), T2
    else
        return ccenergy(T2, ovov)
    end
end

function ccenergy(tiJaB, iajb)
    ecc = 0.0
    nocc = size(tiJaB, 1)
    nvir = size(tiJaB, 4)
    rocc = collect(UnitRange(1, nocc))
    rvir = collect(UnitRange(1, nvir))
    for i in rocc
        for j in rocc
            cache = iajb[i,:,j,:]
            for a in rvir
                for b in rvir
                    ecc += cache[a,b] * 2 * tiJaB[i, j, a, b]
                    ecc -= cache[a,b] * tiJaB[j, i, a, b]
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
                cache = menf[:,e,:,f]
                for n in rocc
                    @simd for m in rocc
                        Fae[a, e] -=
                            cache[m,n] * (2 * tiJaB[m, n, a, f] - tiJaB[n, m, a, f])
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
            cache = menf[:,e,:,f]
            for n in rocc
                for i in rocc
                    @simd for m in rocc
                        Fmi[m, i] +=
                            cache[m,n] * (2 * tiJaB[i, n, e, f] - tiJaB[i, n, f, e])
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
            cache_iajb = iajb[:,a,:,b]
            for j in rocc
                for i in rocc
                    temp = cache_iajb[i, j]
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
    for m in rocc
        for n in rocc
            cache_minj = minj[m,:,n,:]
            cache_menf = menf[m,:,n,:]
            for i in rocc
                for j in rocc
                    #m i n j
                    Wmnij[m, n, i, j] = cache_minj[i, j]
                    for f in rvir
                        for e in rvir
                            Wmnij[m, n, i, j] += tiJaB[i, j, e, f] * cache_menf[e, f] / 2.0
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
            cache_aebf = aebf[:,e,:,f]
            cache_menf = menf[:,e,:,f]
            for b in rvir
                for a in rvir
                    Wabef[a, b, e, f] = cache_aebf[a, b]
                    for n in rocc
                        for m in rocc
                            Wabef[a, b, e, f] += tiJaB[m, n, a, b] * cache_menf[m, n] / 2.0
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
        for m in rocc
            cache_mebj = mebj[m,e,:,:]
            cache_iajb1 = iajb[m,e,:,:]
            cache_iajb2 = iajb[:,e,m,:]
            for b in rvir
                for j in rocc
                    WmBeJ[m, b, e, j] = cache_mebj[b, j]
                    for f in rvir
                        for n in rocc
                            WmBeJ[m, b, e, j] +=
                                cache_iajb1[n, f] *
                                (2 * tiJaB[n, j, f, b] - tiJaB[j, n, f, b]) / 2.0
                            WmBeJ[m, b, e, j] -= cache_iajb2[n, f] * tiJaB[n, j, f, b] / 2.0
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
            cache_mjbe = mjbe[m,:,:,e]
            cache_nemf = nemf[:,e,m,:]
            for b in rvir
                for j in rocc
                    WmBEj[m, b, e, j] = -cache_mjbe[j, b]
                    for f in rvir
                        for n in rocc
                            WmBEj[m, b, e, j] += tiJaB[j, n, f, b] * cache_nemf[n, f] / 2.0
                        end
                    end
                end
            end
        end
    end
    return WmBEj
end
