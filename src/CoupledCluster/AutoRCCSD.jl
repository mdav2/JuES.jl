"""
    JuES.CoupledCluster.AutoRCCSD

Module that performs Restricted CCSD using auto factorized equations.

**Functions:**

    update_energy   Compute CC energy from amplitudes arrays.
    update_amp      Update produces new set of amplitudes from old ones. 
    do_rccsd        Compute RCCSD.

"""
module AutoRCCSD
using JuES
using JuES.Wavefunction
using JuES.IntegralTransformation
using JuES.Output
using TensorOperations
using LinearAlgebra
using Printf

"""
    JuES.CoupledCluster.AutoRCCSD.update_energy(T1::Array{Float64, 2}, T2::Array{Float64, 4}, f::Array{Float64,2}, Voovv::Array{Float64, 4})

Compute CC energy from amplitudes arrays.

**Arguments**

    T1     T1 CC amplitudes array.
    T2     T2 CC amplitudes array.
    f      Fock matrix f(i,a).
    Voovv  ERI tensor V(i,j,a,b).
"""
function update_energy(T1::Array{Float64, 2}, T2::Array{Float64, 4}, f::Array{Float64,2}, Voovv::Array{Float64, 4})

    E::Float64 = 0
    @tensoropt begin
        CC_energy = 2.0*f[k,c]*T1[k,c]
        B[l,c,k,d] := -1.0*T1[l,c]*T1[k,d]
        B[l,c,k,d] += -1.0*T2[l,k,c,d]
        B[l,c,k,d] += 2.0*T2[k,l,c,d]
        CC_energy += B[l,c,k,d]*Voovv[k,l,c,d]
        CC_energy += 2.0*T1[l,c]*T1[k,d]*Voovv[l,k,c,d]
    end
    
    return CC_energy
end

"""
    JuES.CoupledCluster.AutoRCCSD.update_amp(T1::Array{Float64, 2}, T2::Array{Float64, 4}, f::Tuple, V::Tuple, d::Array{Float64, 2}, D::Array{Float64, 4})

Compute new set of T1 and T2 amplitudes from old ones.

**Arguments**

    T1     T1 CC amplitudes array.
    T2     T2 CC amplitudes array.
    f      Tuple containing slices of the full Fock matrix.
    V      Tuple containing slices of the full ERI tensor.
    d      Resolvent tensor D(i,a)
    D      Resolvent tensor D(i,j,a,b)
"""
function update_amp(T1::Array{Float64, 2}, T2::Array{Float64, 4}, f::Tuple, V::Tuple, d::Array{Float64, 2}, D::Array{Float64, 4})

    Voooo, Vooov, Voovv, Vovov, Vovvv, Vvvvv = V
    fock_OO, fock_OV, fock_VV = f

    newT1 = zeros(size(T1))
    newT2 = zeros(size(T2))

    # Get new amplitudes
    @tensoropt (i=>x, j=>x, k=>x, l=>x, a=>100x, b=>100x, c=>100x, d=>100x) begin
        newT1[i,a] += fock_OV[i,a]
        newT1[i,a] -= fock_OO[i,k]*T1[k,a]
        newT1[i,a] += fock_VV[c,a]*T1[i,c]
        newT1[i,a] -= fock_OV[k,c]*T1[i,c]*T1[k,a]
        newT1[i,a] += 2.0*fock_OV[k,c]*T2[i,k,a,c]
        newT1[i,a] -= fock_OV[k,c]*T2[k,i,a,c]
        newT1[i,a] -= T1[k,c]*Vovov[i,c,k,a]
        newT1[i,a] += 2.0*T1[k,c]*Voovv[k,i,c,a]
        newT1[i,a] -= T2[k,i,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] += 2.0*T2[i,k,c,d]*Vovvv[k,a,d,c]
        newT1[i,a] += -2.0*T2[k,l,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += T2[l,k,a,c]*Vooov[k,l,i,c]
        newT1[i,a] += -2.0*T1[k,c]*T1[l,a]*Vooov[l,k,i,c]
        newT1[i,a] -= T1[k,c]*T1[i,d]*Vovvv[k,a,d,c]
        newT1[i,a] += 2.0*T1[k,c]*T1[i,d]*Vovvv[k,a,c,d]
        newT1[i,a] += T1[k,c]*T1[l,a]*Vooov[k,l,i,c]
        newT1[i,a] += -2.0*T1[k,c]*T2[i,l,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[k,c]*T2[l,i,a,d]*Voovv[k,l,c,d]
        newT1[i,a] += T1[k,c]*T2[l,i,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[i,c]*T2[l,k,a,d]*Voovv[l,k,c,d]
        newT1[i,a] += T1[i,c]*T2[l,k,a,d]*Voovv[k,l,c,d]
        newT1[i,a] += -2.0*T1[l,a]*T2[i,k,d,c]*Voovv[k,l,c,d]
        newT1[i,a] += T1[l,a]*T2[i,k,c,d]*Voovv[k,l,c,d]
        newT1[i,a] += T1[k,c]*T1[i,d]*T1[l,a]*Voovv[l,k,c,d]
        newT1[i,a] += -2.0*T1[k,c]*T1[i,d]*T1[l,a]*Voovv[k,l,c,d]
        newT1[i,a] += 4.0*T1[k,c]*T2[i,l,a,d]*Voovv[k,l,c,d]

        newT2[i,j,a,b] += Voovv[i,j,a,b]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*Vvvvv[c,d,a,b]
        newT2[i,j,a,b] += T2[i,j,c,d]*Vvvvv[c,d,a,b]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] += T2[k,l,a,b]*Voooo[i,j,k,l]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,a]*Vovvv[k,b,c,d]
        newT2[i,j,a,b] -= T1[i,c]*T1[j,d]*T1[k,b]*Vovvv[k,a,d,c]
        newT2[i,j,a,b] += T1[i,c]*T1[k,a]*T1[l,b]*Vooov[l,k,j,c]
        newT2[i,j,a,b] += T1[j,c]*T1[k,a]*T1[l,b]*Vooov[k,l,i,c]
        newT2[i,j,a,b] += T2[k,l,a,c]*T2[i,j,d,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[l,j,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[l,k,a,c]*T2[i,j,d,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,d,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T2[i,k,a,c]*T2[l,j,b,d]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,k,a,c]*T2[j,l,b,d]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T2[k,i,a,c]*T2[l,j,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[k,i,a,c]*T2[j,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[i,j,a,c]*T2[l,k,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += -2.0*T2[i,j,a,c]*T2[k,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[k,j,a,c]*T2[i,l,d,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += 4.0*T2[i,k,a,c]*T2[j,l,b,d]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T2[i,j,d,c]*T2[l,k,a,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T1[k,a]*T1[l,b]*Voovv[k,l,c,d]
        newT2[i,j,a,b] += T1[i,c]*T1[j,d]*T2[l,k,a,b]*Voovv[l,k,c,d]
        newT2[i,j,a,b] += T1[k,a]*T1[l,b]*T2[i,j,d,c]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] := -1.0*fock_OO[i,k]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += fock_VV[c,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[k,b]*Vooov[j,i,k,a]
        P_OoVv[i,j,a,b] += T1[j,c]*Vovvv[i,c,a,b]
        P_OoVv[i,j,a,b] += -1.0*fock_OV[k,c]*T1[i,c]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += -1.0*fock_OV[k,c]*T1[k,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,i,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,a]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,b]*Vovov[j,c,k,a]
        P_OoVv[i,j,a,b] += 2.0*T2[i,k,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[i,k,a,c]*Vovov[j,c,k,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,j,a,c]*Vovov[i,c,k,b]
        P_OoVv[i,j,a,b] += -2.0*T1[l,b]*T2[i,k,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[k,i,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,d,b]*Vovvv[k,a,c,d]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[k,i,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,a,d]*Vovvv[k,b,c,d]
        P_OoVv[i,j,a,b] += T1[j,c]*T2[l,k,a,b]*Vooov[l,k,i,c]
        P_OoVv[i,j,a,b] += T1[l,b]*T2[i,k,a,c]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,a]*T2[i,j,d,c]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += T1[k,a]*T2[i,l,c,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += 2.0*T1[j,c]*T2[i,k,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,c]*T2[i,j,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += 2.0*T1[k,c]*T2[i,j,a,d]*Vovvv[k,b,c,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T2[i,l,a,b]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T2[i,l,a,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += T2[j,k,c,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[i,c]*T1[k,a]*T2[j,l,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[k,a]*T2[l,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += T1[i,c]*T1[l,b]*T2[k,j,a,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T2[i,k,d,c]*T2[l,j,a,b]*Voovv[k,l,c,d]
        
        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end

    # Apply the resolvent
    newT1 = newT1.*d
    newT2 = newT2.*D

    # Compute residues 
    r1 = sqrt(sum((newT1 - T1).^2))/length(T1)
    r2 = sqrt(sum((newT2 - T2).^2))/length(T2)

    return newT1, newT2, r1, r2
end

"""
    JuES.CoupledCluster.AutoRCCSD.do_rccsd(wfn::Wfn; kwargs...)

Perform the RCCSD computation.

**Arguments**

    wfn     Wavefunction object.

**Kwargs**

    kwargs...   Options from JuES.
"""
function do_rccsd(wfn::Wfn; kwargs...)

    JuES.CoupledCluster.print_header()
    
    # Process options
    for arg in keys(JuES.CoupledCluster.defaults)
        if arg in keys(kwargs)
            @eval $arg = $(kwargs[arg])
        else
            @eval $arg = $(JuES.CoupledCluster.defaults[arg])
        end
    end

    # Check if the number of electrons is even
    nelec = wfn.nalpha + wfn.nbeta
    nelec % 2 == 0 ? nothing : error("Number of electrons must be even for RHF. Given $nelec")
    nmo = wfn.nmo
    ndocc = Int(nelec/2)
    nvir = nmo - ndocc
    
    # Slices
    o = 1:ndocc
    v = ndocc+1:nmo

    # Get fock matrix
    f = get_fock(wfn, "alpha")

    # Save diagonal terms
    fock_Od = diag(f)[o]
    fock_Vd = diag(f)[v]
    fd = (fock_Od, fock_Vd)

    # Erase diagonal elements from original matrix
    f = f - Diagonal(f)

    # Save useful slices
    fock_OO = f[o,o]
    fock_VV = f[v,v]
    fock_OV = f[o,v]
    f = (fock_OO, fock_OV, fock_VV)

    # Get Necessary ERIs
    V = (get_eri(wfn, "OOOO"), get_eri(wfn, "OOOV"), get_eri(wfn, "OOVV"), get_eri(wfn, "OVOV"), get_eri(wfn, "OVVV"), get_eri(wfn, "VVVV"))
    
    # Auxiliar D matrix
    fock_Od, fock_Vd = fd
    d = [i - a for i = fock_Od, a = fock_Vd]
    d = inv.(d)
    
    D = [i + j - a - b for i = fock_Od, j = fock_Od, a = fock_Vd, b = fock_Vd]
    D = inv.(D)
    
    # Initial Amplitude. Note that f[2] = fock_OV and V[3] = Voovv.
    T1 = f[2].*d
    T2 = D.*V[3]
    
    # Get MP2 energy. Note that f[2] = fock_OV and V[3] = Voovv.
    Ecc = update_energy(T1, T2, f[2], V[3])
    
    @output "MP2 Energy:   {:15.10f}\n\n" Ecc
    
    r1 = 1
    r2 = 1
    dE = 1
    rms = 1
    
    # Start CC iterations
    println("======================================")
    ite = 1
    while abs(dE) > cc_e_conv || rms > cc_max_rms
        if ite > cc_max_iter
            @printf("CC Equations did not converge in %1.0d iterations.", cc_max_iter)
            break
        end
        t = @elapsed begin
            T1, T2, r1, r2 = update_amp(T1, T2, f, V, d, D)
        end
        rms = max(r1,r2)
        oldE = Ecc
        Ecc = update_energy(T1, T2, f[2], V[3])
        dE = Ecc - oldE
        @printf("Iteration %.0f\n", ite)
        @printf("CC Correlation energy: %15.10f\n", Ecc)
        @printf("Energy change:         %15.10f\n", dE)
        @printf("Max RMS residue:       %15.10f\n", rms)
        @printf("Time required:         %15.10f\n", t)
        println("======================================")
        ite += 1
    end
    
    @printf("Final CCSD Energy:     %15.10f\n", Ecc)
end

end #End Module
