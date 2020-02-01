function do_df_rmp2(refWfn::Wfn)
    #build DF-basis
    nocc    = refWfn.nalpha
    nvir    = refWfn.nvira
    C       = refWfn.Ca
    bname   = refWfn.basis.blend()
    df      = psi4.core.BasisSet.build(refWfn.basis.molecule(), "DF_BASIS_MP2", "def2-svp-jkfit")
    dfmints = refWfn.mints
    null    = psi4.core.BasisSet.zero_ao_basis_set()
    pqP     = dfmints.ao_eri(refWfn.basis,refWfn.basis,df,null).np
    Jpq     = dfmints.ao_eri(df,null,df,null).np
    Jpq     = squeeze(Jpq)
    pqP     = squeeze(pqP)
    @inbounds @fastmath begin
        Jpqh    = Jpq^(-1/2)
        _Co      = C[:,1:nocc]
        _Cv    = C[:,nocc+1:nocc+nvir]
    end #@inbounds @fastmath
    @tensoropt begin
        bμν[p,q,Q] := pqP[p,q,P]*Jpqh[P,Q]
    end
    #bμν     = squeeze(bμν)
    @tensoropt begin
        biν[i,ν,Q] := _Co[μ,i]*bμν[μ,ν,Q]
    end
    @tensoropt begin
        bia[i,a,Q] := _Cv[ν,a]*biν[i,ν,Q]
    end
    dmp2 = 0.0
    eps = refWfn.epsa
    bi = zeros(size(bia[1,:,:]))
    bj = zeros(size(bi))
    bAB = zeros(nvir,nvir)
    @fastmath @inbounds for i in 1:nocc
        for j in 1:nocc
            bi[:,:] = bia[i,:,:]
            bj[:,:] = bia[j,:,:]
            @tensoropt begin
                bAB[a,b] = bi[a,Q]*bj[b,Q]
            end
            bBA = transpose(bAB)
            for b in 1:nvir
                for a in 1:nvir
                    iajb = bAB[a,b]
                    ibja = bBA[a,b]
                    dmp2 += iajb*(2*iajb - ibja)/(eps[i] + eps[j] - eps[a+nocc] - eps[b+nocc])
                end
            end
        end
    end
    return dmp2
end

function squeeze(A::AbstractArray)
    #singleton_dims = tuple((d for d in 1:ndims(A) if size(A, d) == 1)...)
    return dropdims(A, dims = (findall(size(A) .== 1)...,))
end
