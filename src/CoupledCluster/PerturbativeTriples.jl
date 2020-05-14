module PerturbativeTriples
using JuES
using JuES.Output
using JuES.Wavefunction
using JuES.IntegralTransformation
using JuES.CoupledCluster: AutoRCCSD
using TensorOperations
using LinearAlgebra

export compute_pT
export do_pT

function do_pT(wfn::Wfn; kwargs...)
    
    Vvvvo, Vvooo, Vvovo = get_eri(wfn, "VVVO"; notation="chem"), get_eri(wfn, "VOOO"; notation="chem"), get_eri(wfn, "VOVO"; notation="chem")

    f = get_fock(wfn, spin="alpha")

    # Save diagonal terms
    fo = diag(f)[1:wfn.nalpha]
    fv = diag(f)[wfn.nalpha+1:wfn.nmo]

    Ecc, T1, T2 = AutoRCCSD.do_rccsd(wfn; kwargs...)
    t = @elapsed begin
        Et = compute_pT(T1, T2, Vvvvo, Vvooo, Vvovo, fo, fv, method="chonky")
    end

    @output "Chonky method:\n"
    @output "CCSD(T) Energy: {:15.10f}\n" Ecc+ Et+wfn.energy
    @output "Time: {:5.5f}\n" t

    t = @elapsed begin
        Et = compute_pT(T1, T2, Vvvvo, Vvooo, Vvovo, fo, fv, method="ijk")
    end

    @output "\nijk method:\n"
    @output "CCSD(T) Energy: {:15.10f}\n" Ecc+ Et+wfn.energy
    @output "Time: {:5.5f}\n" t
end

function compute_pT(T1::Array{Float64, 2}, T2::Array{Float64, 4}, Vvvvo::Array{Float64,4}, Vvooo::Array{Float64,4}, Vvovo::Array{Float64,4}, fo::Array{Float64,1}, fv::Array{Float64,1}; method = "chonky")

    if method == "chonky"
        Et = 0.0
        # Build full resolvent
        Dd = [i + j + k - a - b - c for i = fo, j = fo, k = fo, a = fv, b = fv, c = fv]
        Dd = inv.(Dd)

        @tensoropt (i=>x, j=>x, k=>x, a=>100x, b=>100x, c=>100x) begin
            X[i,j,k,a,b,c] := Vvvvo[b,d,a,i]*T2[k,j,c,d] - Vvooo[c,k,j,l]*T2[i,l,a,b]
            W[i,j,k,a,b,c] := X[i,j,k,a,b,c] + X[i,k,j,a,c,b] + X[k,i,j,c,a,b] + X[k,j,i,c,b,a] + X[j,k,i,b,c,a] + X[j,i,k,b,a,c]
            V[i,j,k,a,b,c] := W[i,j,k,a,b,c] + Vvovo[b,j,c,k]*T1[i,a] + Vvovo[a,i,c,k]*T1[j,b] + Vvovo[a,i,b,j]*T1[k,c]
            X[i,j,k,a,b,c]  = 4*W[i,j,k,a,b,c] + W[i,j,k,b,c,a] + W[i,j,k,c,a,b]
            Y[i,j,k,a,b,c] := V[i,j,k,a,b,c] - V[i,j,k,c,b,a]
        end
        Et = sum(X.*Y.*Dd)/3.0

    elseif method == "ijk"    

        o,v = size(T1)
        Et = 0.0
        for i in 1:o
            for j in 1:o
                for k in 1:o
                    
                    Dd = inv.([fo[i] + fo[j] + fo[k] - a - b - c for a = fv, b = fv, c = fv])

                    # Create views for W
                    Vvvvo_4i = view(Vvvvo, :,:,:,i)
                    Vvvvo_4j = view(Vvvvo, :,:,:,j)
                    Vvvvo_4k = view(Vvvvo, :,:,:,k)

                    Vvooo_2k_3j = view(Vvooo,:,k,j,:)
                    Vvooo_2k_3i = view(Vvooo,:,k,i,:)
                    Vvooo_2i_3j = view(Vvooo,:,i,j,:)
                    Vvooo_2i_3k = view(Vvooo,:,i,k,:)
                    Vvooo_2j_3i = view(Vvooo,:,j,i,:)
                    Vvooo_2j_3k = view(Vvooo,:,j,k,:)

                    T2_1k_2j = view(T2, k,j,:,:)
                    T2_1k_2i = view(T2, k,i,:,:)
                    T2_1i_2j = view(T2, i,j,:,:)
                    T2_1i_2k = view(T2, i,k,:,:)
                    T2_1j_2i = view(T2, j,i,:,:)
                    T2_1j_2k = view(T2, j,k,:,:)

                    T2_1i = view(T2, i,:,:,:)
                    T2_1j = view(T2, j,:,:,:)
                    T2_1k = view(T2, k,:,:,:)

                    # Create views for V
                    T1_1i = view(T1, i, :)
                    T1_1j = view(T1, j, :)
                    T1_1k = view(T1, k, :)

                    Vvovo_2j_4k = view(Vvovo, :,j,:,k)
                    Vvovo_2i_4k = view(Vvovo, :,i,:,k)
                    Vvovo_2i_4j = view(Vvovo, :,i,:,j)

                    @tensoropt begin
                        W[a,b,c] := (Vvvvo_4i[b,d,a]*T2_1k_2j[c,d] - Vvooo_2k_3j[c,l]*T2_1i[l,a,b]  # ijk abc
                                  +  Vvvvo_4i[c,d,a]*T2_1j_2k[b,d] - Vvooo_2j_3k[b,l]*T2_1i[l,a,c]  # ikj acb
                                  +  Vvvvo_4k[a,d,c]*T2_1j_2i[b,d] - Vvooo_2j_3i[b,l]*T2_1k[l,c,a]  # kij cab
                                  +  Vvvvo_4k[b,d,c]*T2_1i_2j[a,d] - Vvooo_2i_3j[a,l]*T2_1k[l,c,b]  # kji cba
                                  +  Vvvvo_4j[c,d,b]*T2_1i_2k[a,d] - Vvooo_2i_3k[a,l]*T2_1j[l,b,c]  # jki bca
                                  +  Vvvvo_4j[a,d,b]*T2_1k_2i[c,d] - Vvooo_2k_3i[c,l]*T2_1j[l,b,a]) # jik bac

                        V[a,b,c] := W[a,b,c] + Vvovo_2j_4k[b,c]*T1_1i[a] + Vvovo_2i_4k[a,c]*T1_1j[b] + Vvovo_2i_4j[a,b]*T1_1k[c]
                        X[a,b,c] := 4*W[a,b,c] + W[b,c,a] + W[c,a,b]
                        Y[a,b,c] := V[a,b,c] - V[c,b,a]
                    end
                    Et += sum(X.*Y.*Dd)/3.0
                end
            end
        end
    end
    return Et
end
end #Module


