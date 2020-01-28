function do_df_rmp2(refWfn::Wfn)
    #build DF-basis
    nocc    = refWfn.nalpha
    nvir    = refWfn.nvira
    C       = refWfn.Ca
    bname   = refWfn.basis.blend()
    df      = psi4.core.BasisSet.build(refWfn.basis.molecule(), "DF_BASIS_MP2", bname)
    #dfmints = psi4.core.MintsHelper(df)
    dfmints = refWfn.mints
    null    = psi4.core.BasisSet.zero_ao_basis_set()
    pqP     = dfmints.ao_eri(refWfn.basis,refWfn.basis,df,null).np
    Jpq     = dfmints.ao_eri(df,null,df,null).np
    Jpq     = squeeze(Jpq)
    pqP     = squeeze(pqP)
    Jpqh    = Jpq^(-1/2)
    @tensor begin
        bμν[p,q,Q] := pqP[p,q,P]*Jpqh[P,Q]
    end
    #bμν     = squeeze(bμν)
    @tensor begin
        biν[i,ν,Q] := C[μ,i]*bμν[μ,ν,Q]
    end
    @tensor begin
        bia[i,a,Q] := C[ν,a]*biν[i,ν,Q]
    end
    @tensor begin
        iajb[i,a,j,b] := bia[i,a,Q]*bia[j,b,Q]
    end
    dmp2 = 0.0
    eps = refWfn.epsa
    for i in 1:nocc
        for j in 1:nocc
            for a in 1:nvir
                for b in 1:nvir
                    aa = a+nocc
                    bb = b+nocc
                    dmp2 += iajb[i,a,j,b]*(2*iajb[i,a,j,b] - iajb[i,b,j,a])/(eps[i] + eps[j] - eps[aa] - eps[bb])
                end
            end
        end
    end
    println(dmp2)
    #for i in 1:nocc
    #    for j in 1:nocc
    #        @tensoropt begin
    #            bAB[a,b] := 
    #        end
end

function squeeze(A::AbstractArray)
    #singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)
    return dropdims(A, dims = (findall(size(A) .== 1)...,))
end
