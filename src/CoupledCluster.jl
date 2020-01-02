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

using Wavefunction
using Base.Threads
using LinearAlgebra
using Dates
export do_rccd


function do_rccd(refWfn::Wfn)
	#implicit maxit = 40
    return do_rccd(refWfn,40)
end
@fastmath @inbounds function do_rccd(refWfn::Wfn,maxit,print=false)
    #goes through appropriate steps to do RCCD
	set_zero_subnormals(true)
    nocc = refWfn.nalpha
    nvir = refWfn.nvira
    iJaB = permutedims(refWfn.pqrs,[1,3,2,4])
	dtt = eltype(iJaB)
    epsa = refWfn.epsa
	T2 = zeros(dtt,nocc,nocc,nvir,nvir)
    Dijab = form_Dijab(T2,epsa)
    T2_init!(T2,iJaB,Dijab)
    Fae = form_Fae(T2,iJaB)
    Fmi = form_Fmi(T2,iJaB)
    Wmnij = form_Wmnij(iJaB,T2)
	Wabef = form_Wabef(iJaB,T2)
	WmBeJ = form_WmBeJ(iJaB,T2)
	WmBEj = form_WmBEj(iJaB,T2)
    for i in UnitRange(1,maxit) #TODO: implement RMS check
        t0 = Dates.Time(Dates.now())
        T2 = cciter(T2,iJaB,Dijab,Fae,Fmi,Wabef,Wmnij,WmBeJ,WmBEj)
        t1 = Dates.Time(Dates.now())
        #print("T2 formed in ")
        #print(convert(Dates.Millisecond, (t1 - t0)))
        #print("\n")
        #t0 = Dates.Time(Dates.now())
		if print
        	print("@CCD ")
        	print(ccenergy(T2,iJaB))
			print("\n")
		end
        t1 = Dates.Time(Dates.now())
        #print("\n")
        #print("energy computed in ")
        #print(convert(Dates.Millisecond, (t1 - t0)))
        #print("\n")
    end
	return ccenergy(T2,iJaB)
end
function ccenergy(tiJaB,iJaB)
	ecc = 0.0
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
	for i in rocc
		for j in rocc
			for a in rvir
				for b in rvir
                    cache = iJaB_oovv[i,j,a,b]
					ecc += cache*2*tiJaB[i,j,a,b]
					ecc -= cache*tiJaB[j,i,a,b]
				end
			end
		end
	end
	return ecc
end

function cciter(tiJaB_i,iJaB,Dijab,Fae,Fmi,Wabef,Wmnij,WmBeJ,WmBEj)
    form_Fae!(Fae,tiJaB_i,iJaB)
    form_Fmi!(Fmi,tiJaB_i,iJaB)
    form_Wmnij!(Wmnij,iJaB,tiJaB_i)
	form_Wabef!(Wabef,iJaB,tiJaB_i)
	form_WmBeJ!(WmBeJ,iJaB,tiJaB_i)
	WmBEj = form_WmBEj!(WmBEj,iJaB,tiJaB_i)
	tiJaB_d = form_T2(tiJaB_i,Fae,Fmi,WmBeJ,WmBEj,Wabef,Wmnij,iJaB,Dijab)
	return tiJaB_d
end

function T2_init!(tiJaB,iJaB,Dijab)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
    tiJaB .= iJaB_oovv ./ Dijab
end

function form_Fae(tiJaB,iJaB)
	dt = eltype(iJaB)
	nvir = size(tiJaB,4)
    Fae = zeros(dt,nvir,nvir)
    form_Fae!(Fae,tiJaB,iJaB)
    return Fae
end
function form_Fae!(Fae,tiJaB,iJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
    Fae .= 0.0
    cache1 = zeros(eltype(iJaB),nocc,nocc)
    cache2 = zeros(eltype(iJaB),nocc,nocc)
    for f in rvir
    	for a in rvir
        	for e in rvir
                for n in rocc
                    @simd for m in rocc
                        Fae[a,e] -= tiJaB[m,n,a,f]*(2*iJaB_oovv[m,n,e,f] - iJaB_oovv[n,m,e,f])
                    end
                end
            end
        end
    end
end

function form_Fmi(tiJaB,iJaB)
	dt = eltype(iJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
    Fmi = zeros(dt,nocc,nocc)
    form_Fmi!(Fmi,tiJaB,iJaB)
	return Fmi
end
function form_Fmi!(Fmi,tiJaB,iJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
    Fmi .= 0.0
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
    for f in rvir
        for e in rvir
            for n in rocc
                for i in rocc
                    @simd for m in rocc
                        Fmi[m,i] += tiJaB[i,n,e,f]*(2*iJaB_oovv[m,n,e,f] - iJaB_oovv[m,n,f,e])
                    end
                end
            end
        end
    end
    #return Fmi
end

function form_Dijab(tiJaB,F)
	dt = eltype(tiJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
	Dijab = zeros(dt,nocc,nocc,nvir,nvir)
	for i in rocc
		for j in rocc
			for a in rvir
				for b in rvir
					aa = a + nocc
					bb = b + nocc
					Dijab[i,j,a,b] = F[i] + F[j] - F[aa] - F[bb]
				end
			end
		end
	end
	return Dijab
end

function form_T2(tiJaB_i,Fae,Fmi,WmBeJ,WmBEj,Wabef,Wmnij,iJaB,Dijab)
	dtt = eltype(tiJaB_i)
	nocc = size(Wmnij,1)
	nvir = size(tiJaB_i,4)
	tiJaB_d = zeros(dtt,nocc,nocc,nvir,nvir)
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
	ttiJaB_i = permutedims(tiJaB_i,[4,1,2,3])
	WWmBEj  = permutedims(WmBEj,[3,1,2,4])
	WWmBeJ  = permutedims(WmBeJ,[3,1,2,4])
	WmBEj = nothing
	WmBeJ = nothing
	Threads.@threads for b in rvir
        for a in rvir 
            _Wabef = Wabef[a,b,:,:]
		   for j in rocc
                for i in rocc
                    temp = iJaB_oovv[i,j,a,b]
                    for e in rvir
                        temp += tiJaB_i[i,j,a,e] * Fae[b,e]
                        temp += tiJaB_i[i,j,e,b] * Fae[a,e]
                        for f in rvir
                            temp += tiJaB_i[i,j,e,f]*_Wabef[e,f]
                        end
                        for m in rocc
                             temp += ttiJaB_i[e,m,i,b]*WWmBEj[e,m,a,j]
                             temp += ttiJaB_i[e,m,j,a]*WWmBEj[e,m,b,i]
                             temp += ttiJaB_i[e,i,m,a]*WWmBeJ[e,m,b,j]
                             temp -= ttiJaB_i[e,m,i,a]*WWmBeJ[e,m,b,j]
                             temp += ttiJaB_i[e,i,m,a]*WWmBeJ[e,m,b,j]
                             temp += ttiJaB_i[e,i,m,a]*WWmBEj[e,m,b,j]
                             temp += ttiJaB_i[e,j,m,b]*WWmBeJ[e,m,a,i]
                             temp -= ttiJaB_i[e,m,j,b]*WWmBeJ[e,m,a,i]
                             temp += ttiJaB_i[e,j,m,b]*WWmBeJ[e,m,a,i]
                             temp += ttiJaB_i[e,j,m,b]*WWmBEj[e,m,a,i]
                        end
                    end

					for m in rocc
                        temp -= tiJaB_i[i,m,a,b]*Fmi[m,j]
                        temp -= tiJaB_i[m,j,a,b]*Fmi[m,i]
                        for n in rocc
                            temp += tiJaB_i[m,n,a,b]*Wmnij[m,n,i,j]
                        end
                    end
                    tiJaB_d[i,j,a,b] += temp
				end
			end
		end
	end
	tiJaB_d .= tiJaB_d ./ Dijab
	return tiJaB_d
end
function form_Wmnij(iJaB,tiJaB)
	dtt = eltype(iJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	Wmnij = zeros(dtt,nocc,nocc,nocc,nocc)
    form_Wmnij!(Wmnij,iJaB,tiJaB)
    return Wmnij
end
@fastmath @inbounds function form_Wmnij!(Wmnij,iJaB,tiJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))

    _tiJaB =  permutedims(tiJaB,[3,4,1,2])
    _iJaB_oovv = permutedims(iJaB_oovv,[3,4,1,2])
	@views Wmnij .= iJaB[1:nocc,1:nocc,1:nocc,1:nocc]
	Threads.@threads for j in rocc
        for n in rocc
	        for i in rocc
	            for m in rocc
	                for f in rvir
	                    @simd for e in rvir
							Wmnij[m,n,i,j] += _tiJaB[e,f,i,j]*_iJaB_oovv[e,f,m,n]/2.0
						end
					end
				end
			end
		end
	end
end

function form_Wabef(iJaB,tiJaB)
	dt = eltype(iJaB)
	nvir = size(tiJaB,4)
	Wabef = zeros(dt,nvir,nvir,nvir,nvir)
    form_Wabef!(Wabef,iJaB,tiJaB)
    return Wabef
end
function form_Wabef!(Wabef,iJaB,tiJaB)
	dtt = eltype(Wabef)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
    #Wabef .= 0
	Wabef .= iJaB[nocc+1:nvir+nocc,nocc+1:nvir+nocc,nocc+1:nvir+nocc,nocc+1:nvir+nocc]
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
    _iJaB_oovv = zeros(dtt,nocc,nocc)
    for f in rvir
        for e in rvir
            _iJaB_oovv .= iJaB_oovv[:,:,e,f]./2.0
	        Threads.@threads for b in rvir
	            for a in rvir
					for n in rocc
						@simd for m in rocc
                            Wabef[a,b,e,f] += tiJaB[m,n,a,b]*_iJaB_oovv[m,n]#*iJaB_oovv[m,n,e,f]
						end
					end
				end
			end
		end
	end
end

function form_WmBeJ(iJaB,tiJaB)
	dtt = eltype(iJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	WmBeJ = zeros(dtt,nocc,nvir,nvir,nocc)
    form_WmBeJ!(WmBeJ,iJaB,tiJaB)
    return WmBeJ
end
function form_WmBeJ!(WmBeJ, iJaB, tiJaB)
	dtt = eltype(WmBeJ)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
    #WmBeJ .= 0 
	@views WmBeJ .= iJaB[1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir,1:nocc]
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
    _iJaB_oovv = zeros(dtt,nocc,nocc)
    for e in rvir
	    for f in rvir
            _iJaB_oovv .= iJaB_oovv[:,:,e,f]./2.0
	        for b in rvir
				for j in rocc
	        		for m in rocc
						@simd for n in rocc
                            WmBeJ[m,b,e,j] -= tiJaB[j,n,f,b]*_iJaB_oovv[m,n]
							WmBeJ[m,b,e,j] += tiJaB[n,j,f,b]*_iJaB_oovv[m,n]*2.0
							WmBeJ[m,b,e,j] -= tiJaB[n,j,f,b]*_iJaB_oovv[n,m]
						end
					end
				end
			end
		end
	end
end


function form_WmBEj(iJaB,tiJaB)
	dtt = eltype(iJaB)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	WmBEj = zeros(dtt,nocc,nvir,nvir,nocc)
    form_WmBEj!(WmBEj,iJaB,tiJaB)
    return WmBEj
end
function form_WmBEj!(WmBEj,iJaB,tiJaB)
	dtt = eltype(WmBEj)
	nocc = size(tiJaB,1)
	nvir = size(tiJaB,4)
	rocc = collect(UnitRange(1,nocc))
	rvir = collect(UnitRange(1,nvir))
	#doing this double permutation scheme reduces floating point performance,
	#but also reduces memory footprint
	iJaB = permutedims(iJaB,[2,1,3,4])
	@views iJaB_ovvo = iJaB[1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir,1:nocc]
	WmBEj -= iJaB_ovvo
	iJaB = permutedims(iJaB,[2,1,3,4])
	@views iJaB_oovv = iJaB[1:nocc,1:nocc,nocc+1:nocc+nvir,nocc+1:nocc+nvir]
    _iJaB_oovv = zeros(dtt,nocc,nocc)
	for e in rvir
		for f in rvir
            _iJaB_oovv .= iJaB_oovv[:,:,e,f]./2.0
        	for m in rocc
	            for b in rvir
				    for j in rocc
					    @simd for n in rocc
                            WmBEj[m,b,e,j] += tiJaB[j,n,f,b]*_iJaB_oovv[n,m]
						end
					end
				end
			end
		end
	end
	return WmBEj
end
end #module CC
