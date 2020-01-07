"""
Module for storing and handling reference wavefunctions.
## Structures
---

Wfn : holds info for disk based and in core computations.
DirectWfn : holds info for integral direct computations.


"""
module Wavefunction

using JuES.DiskTensors
using PyCall
const psi4 = PyNULL()
function __init__()
	copy!(psi4,pyimport("psi4"))
end

export Wfn
export DirectWfn

"""
	Wfn
Data structure for storing integrals, MO coefficients, and misc information about
a reference (HF) wavefunction.

## Fields
---
nalpha::Int number of alpha electrons
nbeta::Int number of beta electrons
nvira::Int number of virtual functions of alpha spin
nvirb::Int number of virtual functions of beta spin
nmo::Int number of molecular orbitals
unrestricted::Bool whether or not the alpha and beta spatial extents are required to
be the same.
Ca::Array{T,2} AO->MO coefficients for alpha MO's
Cb::Array{T,2} AO->MO coefficients for beta MO's
hao::Array{T,2} AO basis core hamiltonian (kinetic + potential)
epsa::Array{T,1} orbital eigenvalues for alpha MO's
epsb::Array{T,1} orbtial eigenvalues for beta MO's
uvsr::Union{Array{T,4},DiskFourTensor} AO basis TEI
pqrs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case AAAA)
pQrS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case ABAB)
pQRs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case ABBA)
PQRS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BBBB)
PqRs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BABA)
PqrS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BAAB)

"""
struct Wfn{T}
	nalpha::Int
	nbeta::Int
    nvira::Int
    nvirb::Int
	nmo::Int
	unrestricted::Bool
    Ca::Array{T,2} #AO->MO coefficients
    Cb::Array{T,2} #AO->MO coefficients
    hao::Array{T,2} #Core hamiltonian
    epsa::Array{T,1} #orbital eigenvalues
    epsb::Array{T,1} #orbital eigenvalues
	uvsr::Union{Array{T,4},DiskFourTensor} #AO basis electron repulsion integrals

	pqrs::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    pQrS::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    pQRs::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    PQRS::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    PqRs::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    PqrS::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
end

"""
	DirectWfn{T}
Data structure for storing information about a wavefunction for which integral
direct procedures will be applied. T::Union{Float32,Float64}

## Fields
---
nalpha::Int number of alpha electrons
nbeta::Int number of beta electrons
nvira::Int number of virtual functions of alpha spin
nvirb::Int number of virtual functions of beta spin
nmo::Int number of molecular orbitals
unrestricted::Bool whether or not the alpha and beta spatial extents are required to
be the same.
Ca::Array{T,2} AO->MO coefficients for alpha MO's
Cb::Array{T,2} AO->MO coefficients for beta MO's
hao::Array{T,2} AO basis core hamiltonian (kinetic + potential)
epsa::Array{T,1} orbital eigenvalues for alpha MO's
epsb::Array{T,1} orbtial eigenvalues for beta MO's
basis::PyObject psi4.core.BasisSet object for accessing basis function info.
mints::PyObject psi4.core.MintsHelper object for computing integrals.
"""
struct DirectWfn{T}
	#requires the module using this struct to have properly imported
	#psi4 with pycall
	nalpha::Int
	nbeta::Int
    nvira::Int
    nvirb::Int
	nmo::Int
	unrestricted::Bool
    Ca::Array{T,2} #AO->MO coefficients
    Cb::Array{T,2} #AO->MO coefficients
    hao::Array{T,2} #Core hamiltonian
    epsa::Array{T,1} #orbital eigenvalues
    epsb::Array{T,1} #orbital eigenvalues
	basis::PyObject
	mints::PyObject
end
function Wfn(wfn::PyObject)
	Wfn(wfn,Float64,false,false)
end

function Wfn(wfn,dt,unrestricted::Bool,diskbased::Bool)
    dummy2 = Array{dt}(undef,0,0) #placeholder 2D array
    dummy4 = Array{dt}(undef,0,0,0,0) #placeholder 4D array
    mints = psi4.core.MintsHelper(wfn.basisset()) #to generate integrals
    nbf   = wfn.nmo()
    nocca =  wfn.nalpha()
    nvira = nbf - nocca
    noccb = wfn.nbeta()
    nvirb = nbf - noccb
    epsa  = convert(Array{dt,1},wfn.epsilon_a().to_array()) #orbital eigenvalues
    epsb  = convert(Array{dt,1},wfn.epsilon_b().to_array()) #orbital eigenvalues
    _Ca   = wfn.Ca() #as psi4 Matrix objects for use with MintsHelper
    _Cb   = wfn.Cb() #as psi4 Matrix objects for use with MintsHelper
    Ca    = convert(Array{dt,2}, _Ca.to_array())
    Cb    = convert(Array{dt,2},_Cb.to_array())
    hao   = convert(Array{dt,2}, wfn.H().to_array()) #core hamiltonian in AO
    uvsr  = convert(Array{dt,4},mints.ao_eri().to_array()) #AO basis integrals
    pqrs  = convert(Array{dt,4},mints.mo_eri(_Ca,_Ca,_Ca,_Ca).to_array()) #MO basis integrals
    if unrestricted #avoid making these if not an unrestricted or open shell wfn
		#various spin cases notation --> alpha BETA
        pQrS  = convert(Array{dt,4},mints.mo_eri(_Ca,_Ca,_Cb,_Cb).to_array())
        pQRs  = convert(Array{dt,4},mints.mo_eri(_Ca,_Cb,_Cb,_Ca).to_array())
        PQRS  = convert(Array{dt,4},mints.mo_eri(_Cb,_Cb,_Cb,_Cb).to_array())
        PqRs  = convert(Array{dt,4},mints.mo_eri(_Cb,_Ca,_Cb,_Ca).to_array())
        PqrS  = convert(Array{dt,4},mints.mo_eri(_Cb,_Ca,_Ca,_Cb).to_array())
    else
		#just fill with placeholder for RHF case
        pQrS  = dummy4
        pQRs  = dummy4
        PQRS  = dummy4
        PqRs  = dummy4
        PqrS  = dummy4
    end
	#create the Wfn object and return it!
    owfn = Wfn{dt}(nocca,noccb,nvira,nvirb,nbf,unrestricted,
           Ca,Cb,hao,
           epsa,epsb,
           uvsr,pqrs,pQrS,pQRs,PQRS,PqRs,PqrS)
    return owfn
end
function DirectWfn(wfn)
	"""
	Constructor for DirectWfn objects.
	"""
    _Ca   = wfn.Ca() #as psi4 Matrix objects for use with MintsHelper
    _Cb   = wfn.Cb() #as psi4 Matrix objects for use with MintsHelper
    nbf   = wfn.nmo()
    nalpha=  wfn.nalpha()
    nbeta = wfn.nbeta()
    nvira = nbf - nalpha
    nvirb = nbf - nbeta
	nmo   = wfn.nmo()
	unrestricted = false
	dt = Float64 #TODO: make type flexible for Float32 support
    Ca    = convert(Array{Float64,2}, _Ca.to_array())
    Cb    = convert(Array{dt,2},_Cb.to_array())
    hao   = convert(Array{Float64,2}, wfn.H().to_array()) #core hamiltonian in AO
    epsa  = convert(Array{dt,1},wfn.epsilon_a().to_array()) #orbital eigenvalues
    epsb  = convert(Array{dt,1},wfn.epsilon_b().to_array()) #orbital eigenvalues
	basis = wfn.basisset()
    mints = psi4.core.MintsHelper(wfn.basisset()) #to generate integrals
	DirectWfn(nalpha,nbeta,nvira,nvirb,nmo,unrestricted,Ca,Cb,hao,epsa,epsb,basis,mints)
end
end
