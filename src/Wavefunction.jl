"""
Module for storing and handling reference wavefunctions.

## structs

    Wfn -> holds info for disk based and in core computations.
    DirectWfn -> holds info for integral direct computations.

"""
module Wavefunction

using JuES.DiskTensors
using JuES.Transformation
using JuES.Integrals
using JuES
using PyCall


export Wfn
export DirectWfn

"""
	Wfn
Data structure for storing integrals, MO coefficients, and misc information about
a reference (HF) wavefunction.

## Fields
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
    basis::PyObject
    mints::PyObject
    Ca::Array{T,2} #AO->MO coefficients
    Cb::Array{T,2} #AO->MO coefficients
    Cao::Array{T,2}
    Cav::Array{T,2}
    Cbo::Array{T,2}
    Cbv::Array{T,2}
    hao::Array{T,2} #Core hamiltonian
    epsa::Array{T,1} #orbital eigenvalues
    epsb::Array{T,1} #orbital eigenvalues
    uvsr::Union{Array{T,4},DiskFourTensor} #AO basis electron repulsion integrals

    ijab::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    iJaB::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    iJAb::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    IJAB::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    IjAb::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    IjaB::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
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
    Wfn(wfn, Float64, false, false)
end

function Wfn(wfn, dt, unrestricted::Bool, diskbased::Bool, name::String = "default")
    dummy2 = Array{dt}(undef, 0, 0) #placeholder 2D array
    dummy4 = Array{dt}(undef, 0, 0, 0, 0) #placeholder 4D array
    basis = wfn.basisset()
    mints = psi4.core.MintsHelper(basis) #to generate integrals
    nbf = wfn.nmo()
    nocca = wfn.nalpha()
    nvira = nbf - nocca
    noccb = wfn.nbeta()
    nvirb = nbf - noccb
    epsa = convert(Array{dt,1}, wfn.epsilon_a().to_array()) #orbital eigenvalues
    epsb = convert(Array{dt,1}, wfn.epsilon_b().to_array()) #orbital eigenvalues
    _Ca = wfn.Ca() #as psi4 Matrix objects for use with MintsHelper
    _Cb = wfn.Cb() #as psi4 Matrix objects for use with MintsHelper
    _Cao = wfn.Ca_subset("AO", "OCC")
    _Cav = wfn.Ca_subset("AO", "VIR")
    _Cbo = wfn.Cb_subset("AO", "OCC")
    _Cbv = wfn.Cb_subset("AO", "VIR")
    Cao = convert(Array{dt,2}, wfn.Ca_subset("AO", "OCC").to_array())
    Cav = convert(Array{dt,2}, wfn.Ca_subset("AO", "VIR").to_array())
    Cbo = convert(Array{dt,2}, wfn.Cb_subset("AO", "OCC").to_array())
    Cbv = convert(Array{dt,2}, wfn.Cb_subset("AO", "VIR").to_array())
    Ca = convert(Array{dt,2}, _Ca.to_array())
    Cb = convert(Array{dt,2}, _Cb.to_array())
    hao = convert(Array{dt,2}, wfn.H().to_array()) #core hamiltonian in AO
    if diskbased
        uvsr = disk_ao(mints, basis, string("jues.", "$name", ".uvsr.0"))
    else
        uvsr = convert(Array{dt,4}, mints.ao_eri().to_array()) #AO basis integrals
    end
    if diskbased
        ijab = tei_transform(uvsr, Cao, Cav, Cao, Cav, "ijab")
    else
        #pqrs  = convert(Array{dt,4},mints.mo_eri(_Ca,_Ca,_Ca,_Ca).to_array()) #MO basis integrals
        ijab = tei_transform(uvsr, Cao, Cav, Cao, Cav, "ijab")
    end
    if unrestricted #avoid making these if not an unrestricted or open shell wfn
        #various spin cases notation --> alpha BETA
        if diskbased
            iJaB = tei_transform(uvsr, Cao, Cav, Cbo, Cbv, string("$name", ".iJaB"))
            iJAb = tei_transform(uvsr, Cao, Cbv, Cbo, Cav, string("$name", ".iJAb"))
            IJAB = tei_transform(uvsr, Cbo, Cbv, Cbo, Cbv, string("$name", ".IJAB"))
            IjAb = tei_transform(uvsr, Cbo, Cbv, Cao, Cav, string("$name", ".IjAb"))
            IjaB = tei_transform(uvsr, Cbo, Cav, Cao, Cbv, string("$name", ".IjaB"))
        else
            iJaB = tei_transform(uvsr, Cao, Cav, Cbo, Cbv, string("$name", ".iJaB"))
            iJAb = tei_transform(uvsr, Cao, Cbv, Cbo, Cav, string("$name", ".iJAb"))
            IJAB = tei_transform(uvsr, Cbo, Cbv, Cbo, Cbv, string("$name", ".IJAB"))
            IjAb = tei_transform(uvsr, Cbo, Cbv, Cao, Cav, string("$name", ".IjAb"))
            IjaB = tei_transform(uvsr, Cbo, Cav, Cao, Cbv, string("$name", ".IjaB"))
        end
    else
        #just fill with placeholder for RHF case
        iJaB = dummy4
        iJAb = dummy4
        IJAB = dummy4
        IjAb = dummy4
        IjaB = dummy4
    end
    #create the Wfn object and return it!
    owfn = Wfn{dt}(
        nocca,
        noccb,
        nvira,
        nvirb,
        nbf,
        unrestricted,
        basis,
        mints,
        Ca,
        Cb,
        Cao,
        Cav,
        Cbo,
        Cbv,
        hao,
        epsa,
        epsb,
        uvsr,
        ijab,
        iJaB,
        iJAb,
        IJAB,
        IjAb,
        IjaB,
    )
    return owfn
end
function DirectWfn(wfn)
    """
    Constructor for DirectWfn objects.
    """
    _Ca = wfn.Ca() #as psi4 Matrix objects for use with MintsHelper
    _Cb = wfn.Cb() #as psi4 Matrix objects for use with MintsHelper
    nbf = wfn.nmo()
    nalpha = wfn.nalpha()
    nbeta = wfn.nbeta()
    nvira = nbf - nalpha
    nvirb = nbf - nbeta
    nmo = wfn.nmo()
    unrestricted = false
    dt = Float64 #TODO: make type flexible for Float32 support
    Ca = convert(Array{Float64,2}, _Ca.to_array())
    Cb = convert(Array{dt,2}, _Cb.to_array())
    hao = convert(Array{Float64,2}, wfn.H().to_array()) #core hamiltonian in AO
    epsa = convert(Array{dt,1}, wfn.epsilon_a().to_array()) #orbital eigenvalues
    epsb = convert(Array{dt,1}, wfn.epsilon_b().to_array()) #orbital eigenvalues
    basis = wfn.basisset()
    mints = psi4.core.MintsHelper(wfn.basisset()) #to generate integrals
    DirectWfn(
        nalpha,
        nbeta,
        nvira,
        nvirb,
        nmo,
        unrestricted,
        Ca,
        Cb,
        hao,
        epsa,
        epsb,
        basis,
        mints,
    )
end
end
