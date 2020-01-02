module CISingles
"""
Module for performing CIS computations.

Currently only RHF reference supported
"""

using IterativeSolvers
using JuES.Wavefunction
using PyCall
using LinearAlgebra 
using Dates

const psi4 = PyNULL()

function __init__()
	copy!(psi4,pyimport("psi4"))
end
dt = Float64

export setup_rcis
export do_CIS
@inline @inbounds @fastmath function HS(eri::Array{dt,4}, F::Array{dt,2}, i::Int64, j::Int64, a::Int64, b::Int64)
	"computes a singly excited Hamiltonian matrix element. Slaters rules are
	incorporated via the kronecker deltas"
    return kron(i,j)*F[a,b] - kron(a,b)*F[i,j] + so_eri(eri,a,j,i,b)
end

@inbounds @fastmath function make_H(nocc,nvir,so_eri,F)
	"Constructs Hamiltonian matrix"
    H = zeros(dt,nocc*nvir,nocc*nvir)
    rocc = collect(UnitRange(1,nocc))
    rvir = collect(UnitRange(1,nvir))
    for b in rvir
        for i in rocc
            for j in rocc
                for a in rvir
                    aa = a + nocc
                    bb = b + nocc
                    I = (i-1)*nvir+a
                    J = (j-1)*nvir+b
                    H[I,J] = HS(so_eri,F,i,j,aa,bb)
                end
            end
        end
    end
    return H
end
function setup_rcis(wfn,dt)
    _C = wfn.Ca()
    nbf = wfn.nmo()
    #nso = 2*nbf
    #nocc = 2*wfn.nalpha()
    nso = nbf
    nocc = wfn.nalpha()
    eps = wfn.epsilon_a().to_array()
    nvir = nso - nocc
    basis = wfn.basisset()
    mints = psi4.core.MintsHelper(basis)
    mo_eri = convert(Array{dt},mints.mo_eri(_C,_C,_C,_C).to_array())
    ff = zeros(dt,nso,nso)
    r = collect(UnitRange(1,nso))
    @inbounds @fastmath for i in r
        ff[i,i] = eps[Int64(fld((i+1),2))]
    end
    return nocc,nvir,mo_eri,ff
end


function do_CIS(refWfn::Wfn,nroots,algo="lobpcg",doprint=false)
	println("in do_CIS")
	nocc = 2*refWfn.nalpha
	nvir = 2*refWfn.nvira
	so_eri = refWfn.pqrs
	F = zeros(nocc + nvir, nocc + nvir)
	r = collect(UnitRange(1,nocc+nvir))
    @inbounds @fastmath for i in r
        F[i,i] = refWfn.epsa[Int64(fld((i+1),2))]
    end
    t0 = Dates.Time(Dates.now())
    t1 = Dates.Time(Dates.now())
    if doprint
        println("constructing Hamiltonian ") 
    end
    t0 = Dates.Time(Dates.now())
    H = make_H(nocc,nvir,so_eri,F)
    t1 = Dates.Time(Dates.now())
    if doprint
        print("Hamiltonian constructed in ") 
        print(convert(Dates.Millisecond,(t1 - t0)))
        print("\n")
    end
    if algo == "lobpcg" || algo == "iter"
        t0 = Dates.Time(Dates.now())
        eigs = lobpcg(H,false,nroots).λ
        t1 = Dates.Time(Dates.now())
        if doprint
            print("Hamiltonian iteratively solved in ") 
            print(convert(Dates.Millisecond,(t1 - t0)))
            print("\n")
        end
	elseif algo == "davidson"
		println(size(H))
        t0 = Dates.Time(Dates.now())
        eigs = eigdav(H,1,4,100,1E-6)
        t1 = Dates.Time(Dates.now())
		if doprint
            print("Hamiltonian iteratively solved in ") 
            print(convert(Dates.Millisecond,(t1 - t0)))
            print("\n")
		end
	elseif algo == "svd"
        t0 = Dates.Time(Dates.now())
		eigs = svdvals(H)
		print(eigs)
        t1 = Dates.Time(Dates.now())
        if doprint
            print("Hamiltonian diagonalized exactly in ") 
            print(convert(Dates.Millisecond,(t1 - t0)))
            print("\n")
        end
    elseif algo == "diag"
    	H = Symmetric(H)
        t0 = Dates.Time(Dates.now())
        eigs = eigvals(H,1:nroots)
        t1 = Dates.Time(Dates.now())
        if doprint
            print("Hamiltonian diagonalized exactly in ") 
            print(convert(Dates.Millisecond,(t1 - t0)))
            print("\n")
        end
    else
        if doprint
            println("solver ",algo," is not supported! Choose from { lobpcg, diag }")
        end
        return false
    end
    return eigs[1:nroots]
end

@views @inline @fastmath @inbounds function so_eri(mo_eri,p,r,q,s)
    pp = Int64(fld((p+1),2))
    qq = Int64(fld((q+1),2))
    rr = Int64(fld((r+1),2))
    ss = Int64(fld((s+1),2))
    return (p%2==q%2)*(r%2==s%2)*mo_eri[pp,qq,rr,ss] - (p%2==s%2)*(q%2==r%2)*mo_eri[pp,ss,qq,rr]
end

@inline @fastmath function kron(a::Int64,b::Int64)
    "Kronecker delta function"
    @fastmath return a == b
end
end
