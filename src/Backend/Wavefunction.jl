"""
Module for storing and handling reference wavefunctions.

## structs

    Wfn -> holds info for disk based and in core computations.
    DirectWfn -> holds info for integral direct computations.

"""
module Wavefunction

using JuES.DiskTensors
using JuES.Transformation
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

ao_eri::Union{Array{T,4},DiskFourTensor} AO basis TEI

pqrs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case AAAA)

pQrS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case ABAB)

pQRs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case ABBA)

PQRS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BBBB)

PqRs::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BABA)

PqrS::Union{Array{T,4},DiskFourTensor} MO basis TEI (spin case BAAB)
"""
struct Wfn{T}
    vnuc
    nalpha::Int
    nbeta::Int
    nvira::Int
    nvirb::Int
    nmo::Int
    unrestricted::Bool
    basis
    mints
    Ca::Array{T,2} #AO->MO coefficients
    Cb::Array{T,2} #AO->MO coefficients
    Cao::Array{T,2}
    Cav::Array{T,2}
    Cbo::Array{T,2}
    Cbv::Array{T,2}
    hao::Array{T,2} #Core hamiltonian
    epsa::Array{T,1} #orbital eigenvalues
    epsb::Array{T,1} #orbital eigenvalues
    ao_eri::Union{Array{T,4},DiskFourTensor} #AO basis electron repulsion integrals

    #ijab::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    #iJaB::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    #iJAb::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    #IJAB::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    #IjAb::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
    #IjaB::Union{Array{T,4},DiskFourTensor} #MO basis electron repulsion integrals
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
    Wfn{Float64}(wfn)
end

"""
    Wfn{T}(wfn::PyObject; unrestricted::Bool=false, diskbased::Bool=false, name::String = "default", df::Bool=false) where T
"""
function Wfn{T}(wfn::PyObject; unrestricted::Bool=false, diskbased::Bool=false, name::String = "default", df::Bool=false) where T
    dt = T
    dummy2 = Array{dt}(undef, 0, 0) #placeholder 2D array
    dummy4 = Array{dt}(undef, 0, 0, 0, 0) #placeholder 4D array
    vnuc = wfn.molecule().nuclear_repulsion_energy()
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
        ao_eri = disk_ao(mints, basis, string("jues.", "$name", ".ao_eri.0"))
    elseif !df
        ao_eri = convert(Array{dt,4}, mints.ao_eri().to_array()) #AO basis integrals
    end
    if diskbased
        ijab = tei_transform(ao_eri, Cao, Cav, Cao, Cav, "ijab")
    else
        #pqrs  = convert(Array{dt,4},mints.mo_eri(_Ca,_Ca,_Ca,_Ca).to_array()) #MO basis integrals
        ijab = tei_transform(ao_eri, Cao, Cav, Cao, Cav, "ijab")
    end
    if unrestricted #avoid making these if not an unrestricted or open shell wfn
        #various spin cases notation --> alpha BETA
        if diskbased
            iJaB = tei_transform(ao_eri, Cao, Cav, Cbo, Cbv, string("$name", ".iJaB"))
            iJAb = tei_transform(ao_eri, Cao, Cbv, Cbo, Cav, string("$name", ".iJAb"))
            IJAB = tei_transform(ao_eri, Cbo, Cbv, Cbo, Cbv, string("$name", ".IJAB"))
            IjAb = tei_transform(ao_eri, Cbo, Cbv, Cao, Cav, string("$name", ".IjAb"))
            IjaB = tei_transform(ao_eri, Cbo, Cav, Cao, Cbv, string("$name", ".IjaB"))
        else
            iJaB = tei_transform(ao_eri, Cao, Cav, Cbo, Cbv, string("$name", ".iJaB"))
            iJAb = tei_transform(ao_eri, Cao, Cbv, Cbo, Cav, string("$name", ".iJAb"))
            IJAB = tei_transform(ao_eri, Cbo, Cbv, Cbo, Cbv, string("$name", ".IJAB"))
            IjAb = tei_transform(ao_eri, Cbo, Cbv, Cao, Cav, string("$name", ".IjAb"))
            IjaB = tei_transform(ao_eri, Cbo, Cav, Cao, Cbv, string("$name", ".IjaB"))
        end
    else
        #just fill with placeholder for RHF case
        iJaB = dummy4
        iJAb = dummy4
        IJAB = dummy4
        IjAb = dummy4
        IjaB = dummy4
    end
    # create the Wfn object and return it!
    owfn = Wfn{dt}(
        vnuc,
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
        ao_eri 
        #ijab,
        #iJaB,
        #iJAb,
        #IJAB,
        #IjAb,
        #IjaB,
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

"""
    disk_ao

computes all AO basis TEI and fills a DiskFourTensor object with those
"""
function disk_ao(mints::PyObject, basis::PyObject, name::String = "default")
    integ = mints.integral()
    si = integ.shells_iterator()
    si.first()
    nao = basis.nbf()
    ao_eri = DiskFourTensor("/tmp/disk_gao.$name.jues.0", Float64, nao, nao, nao, nao, "w")
    blockfill!(ao_eri, 0.0)
    while !si.is_done()
        p, q, r, s = (si.p, si.q, si.r, si.s)
        pf = basis.shell_to_basis_function(p) + 1
        qf = basis.shell_to_basis_function(q) + 1
        rf = basis.shell_to_basis_function(r) + 1
        sf = basis.shell_to_basis_function(s) + 1
        shell = mints.ao_eri_shell(p, q, r, s).to_array()
        pn, qn, rn, sn = size(shell)
        pn -= 1
        qn -= 1
        rn -= 1
        sn -= 1
        ao_eri[pf:pf+pn, qf:qf+qn, rf:rf+rn, sf:sf+sn] = shell
        ao_eri[qf:qf+qn, pf:pf+pn, rf:rf+rn, sf:sf+sn] = permutedims(shell, [2, 1, 3, 4])
        ao_eri[pf:pf+pn, qf:qf+qn, sf:sf+sn, rf:rf+rn] = permutedims(shell, [1, 2, 4, 3])
        ao_eri[qf:qf+qn, pf:pf+pn, sf:sf+sn, rf:rf+rn] =
            permutedims(permutedims(shell, [2, 1, 3, 4]), [1, 2, 4, 3])
        ao_eri[rf:rf+rn, sf:sf+sn, pf:pf+pn, qf:qf+qn] =
            permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2])
        ao_eri[sf:sf+sn, rf:rf+rn, pf:pf+pn, qf:qf+qn] = permutedims(
            permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2]),
            [2, 1, 3, 4],
        )
        ao_eri[rf:rf+rn, sf:sf+sn, qf:qf+qn, pf:pf+pn] = permutedims(
            permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2]),
            [1, 2, 4, 3],
        )
        ao_eri[sf:sf+sn, rf:rf+rn, qf:qf+qn, pf:pf+pn] =
            permutedims(permutedims(shell, [4, 2, 3, 1]), [1, 3, 2, 4])
        si.next()
    end
    return ao_eri

end

end # Module
