module UNFRCCSD
using JuES.Wavefunction
using JuES.Transformation
using TensorOperations
using LinearAlgebra
using Printf

function update_energy(T1::Array{Float64, 2}, T2::Array{Float64, 4}, f::Array{Float64,2}, V::Array{Float64, 4})

    E::Float64 = 0
    @tensoropt begin
        CC_energy = 2.0*f[k,c]*T1[k,c]
        B_OVOV[l,c,k,d] := -1.0*T1[l,c]*T1[k,d]
        B_OVOV[l,c,k,d] += -1.0*T2[l,k,c,d]
        B_OVOV[l,c,k,d] += 2.0*T2[k,l,c,d]
        CC_energy += B_OVOV[l,c,k,d]*V[k,l,c,d]
        CC_energy += 2.0*T1[l,c]*T1[k,d]*V[l,k,c,d]
    end
    
    return CC_energy
end

function update_amp(T1::Array{Float64, 2}, T2::Array{Float64, 4}, f::Tuple, V::Tuple, d::Array{Float64, 2}, D::Array{Float64, 4}, info::Dict)

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
        P_OoVv[i,j,a,b] += 1.0*fock_VV[c,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[k,b]*Vooov[j,i,k,a]
        P_OoVv[i,j,a,b] += 1.0*T1[j,c]*Vovvv[i,c,a,b]
        P_OoVv[i,j,a,b] += -1.0*fock_OV[k,c]*T1[i,c]*T2[k,j,a,b]
        P_OoVv[i,j,a,b] += -1.0*fock_OV[k,c]*T1[k,a]*T2[i,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,i,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,a]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T1[i,c]*T1[k,b]*Vovov[j,c,k,a]
        P_OoVv[i,j,a,b] += 2.0*T2[i,k,a,c]*Voovv[k,j,c,b]
        P_OoVv[i,j,a,b] += -1.0*T2[i,k,a,c]*Vovov[j,c,k,b]
        P_OoVv[i,j,a,b] += -1.0*T2[k,j,a,c]*Vovov[i,c,k,b]
        P_OoVv[i,j,a,b] += -2.0*T1[l,b]*T2[i,k,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += 1.0*T1[l,b]*T2[k,i,a,c]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,d,b]*Vovvv[k,a,c,d]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[k,i,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += -1.0*T1[j,c]*T2[i,k,a,d]*Vovvv[k,b,c,d]
        P_OoVv[i,j,a,b] += 1.0*T1[j,c]*T2[l,k,a,b]*Vooov[l,k,i,c]
        P_OoVv[i,j,a,b] += 1.0*T1[l,b]*T2[i,k,a,c]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,a]*T2[i,j,d,c]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += 1.0*T1[k,a]*T2[i,l,c,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += 2.0*T1[j,c]*T2[i,k,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += -1.0*T1[k,c]*T2[i,j,a,d]*Vovvv[k,b,d,c]
        P_OoVv[i,j,a,b] += 2.0*T1[k,c]*T2[i,j,a,d]*Vovvv[k,b,c,d]
        P_OoVv[i,j,a,b] += 1.0*T1[k,c]*T2[i,l,a,b]*Vooov[k,l,j,c]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T2[i,l,a,b]*Vooov[l,k,j,c]
        P_OoVv[i,j,a,b] += 1.0*T2[j,k,c,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += 1.0*T1[k,c]*T1[j,d]*T2[i,l,a,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += 1.0*T1[k,c]*T1[l,a]*T2[i,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += 1.0*T1[i,c]*T1[k,a]*T2[l,j,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T1[i,c]*T1[k,a]*T2[j,l,b,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += 1.0*T1[i,c]*T1[k,a]*T2[l,j,d,b]*Voovv[l,k,c,d]
        P_OoVv[i,j,a,b] += 1.0*T1[i,c]*T1[l,b]*T2[k,j,a,d]*Voovv[k,l,c,d]
        P_OoVv[i,j,a,b] += -2.0*T2[i,k,d,c]*T2[l,j,a,b]*Voovv[k,l,c,d]
        
        newT2[i,j,a,b] += P_OoVv[i,j,a,b] + P_OoVv[j,i,b,a]
    end

    newT1 = newT1.*d
    newT2 = newT2.*D

    r1 = sqrt(sum((newT1 - T1).^2))/length(T1)
    r2 = sqrt(sum((newT2 - T2).^2))/length(T2)

    return newT1, newT2, r1, r2
end


function get_integrals(wfn::Wfn, info::Dict)

    C = wfn.Ca
    gao = wfn.uvsr
    @tensoropt begin
        Vchem[i,a,j,b] := C[σ,b]*C[λ,j]*C[ν,a]*C[μ,i]*gao[μ,ν,λ,σ]
    end
    hao = wfn.hao

    # Slices
    o = 1:info["ndocc"]
    v = info["ndocc"]+1:info["nmo"]

    # Auxiliar slices for the F matrix
    A = view(Vchem, :, :, o, o)
    B = view(Vchem, :, o, :, o)

    # Form the full Fock matrices
    @tensoropt begin
        h[p,q]  := C[u,p]*C[v,q]*hao[u,v]
        Va[p,q] := A[p,q,k,k]
        Vb[p,q] := B[p,k,q,k]
    end

    f  = h + 2Va - Vb

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

    # Save slices of two-electron repulsion integral
    @tensor begin
        Vphys[p,r,q,s] := Vchem[p,q,r,s]
    end

    V = (Vphys[o,o,o,o], Vphys[o,o,o,v], Vphys[o,o,v,v], Vphys[o,v,o,v], Vphys[o,v,v,v], Vphys[v,v,v,v])

    return fd, f, V

end

function do_rccsd(wfn::Wfn)
info = Dict()
# Check if the number foe electrons is even
info["nelec"] = wfn.nalpha + wfn.nbeta
info["nelec"] % 2 == 0 ? nothing : throw(DomainError(info["nelec"], "Number of electrons must be even for RHF"))
info["nmo"] = wfn.nmo
info["ndocc"] = Int(info["nelec"]/2)
info["nvir"] = info["nmo"] - info["ndocc"]


println("Number of electrons:              $(info["nelec"])")
println("Number of Doubly Occupied MOs:    $(info["ndocc"])")
println("Number of MOs:                    $(info["nmo"])")

fd, f, V = get_integrals(wfn, info)

# Auxiliar D matrix
fock_Od, fock_Vd = fd
d = [i - a for i = fock_Od, a = fock_Vd]
d = inv.(d)

D = [i + j - a - b for i = fock_Od, j = fock_Od, a = fock_Vd, b = fock_Vd]
D = inv.(D)

# Initial Amplitude
T1 = f[2].*d
T2 = D.*V[3]

# Get MP2 energy

Ecc = update_energy(T1, T2, f[2], V[3])

@printf("MP2 Energy:   %15.10f\n", Ecc)

# Set up iteration options
r1 = 1
r2 = 1
dE = 1
ite = 1
rms_LIM = 10^-8
E_LIM = 10^-12
rms = 1

println("======================================")


# Start CC iterations

while abs(dE) > E_LIM || rms > rms_LIM
    if ite > 50
        break
    end
    t = @elapsed begin
        T1, T2, r1, r2 = update_amp(T1, T2, f, V, d, D, info)
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

@printf("CCSD Energy:     %15.10f\n", Ecc)
end
end #End Module
