function do_direct_rmp2(refWfn::Wfn)
    dmp2 = 0.0
    nocc = refWfn.nalpha
    rocc = 1:1:refWfn.nalpha
    epsa = refWfn.epsa
    Ca = refWfn.Ca
    mints = refWfn.mints
    basis = refWfn.basis
    nao = basis.nbf()
    nvir = nao - nocc
    rvir = 1:nvir
    eps = refWfn.epsa
    dmp2 = 0.0
    for i in rocc
        Cav = Ca[:,nocc+1:nao]
        Cao = Ca[:,1:nocc]
        D = zeros(nvir,nocc,nvir)
        #out = zeros(size(D))
        for j in rocc
            for a in rvir
                for b in rvir
                    aa = a + nocc
                    bb = b + nocc
                    D[a,j,b] = 1/(eps[i] + eps[j] - eps[aa] - eps[bb])
                end
            end
        end
        #println("denominator formed in $time")
        temp = direct_ao_contract(i,Ca,mints,basis)
        #println("integrals computed in $time")
        @tensoropt begin
            temp2[a,λ,σ] := Cav[ν,a]*temp[ν,λ,σ]
        end
        temp = nothing
        @tensoropt begin
            temp3[a,j,σ] := Cao[λ,j]*temp2[a,λ,σ]
        end
        temp2 = nothing
        @tensoropt begin
            temp4[a,j,b] := Cav[σ,b]*temp3[a,j,σ]
        end
        temp3 = nothing
        out = temp4 .* (2*temp4 - permutedims(temp4,[3,2,1])) .* D
        #@tensor begin
        #    out[a,j,b] := temp4[a,j,b] * (2*temp4[a,j,b] - temp4[b,j,a]) * D[a,j,b]
        #end
        dmp2 += reduce(+,out)
    end
    return dmp2

end
"""
contraction to N^3 (iν|λσ) given specific i
"""
function direct_ao_contract(i,C::Array{Float64,2},mints::PyObject, basis::PyObject)
    #integ = mints.integral()
    #si = integ.shells_iterator()
    #si.first()
    #eri = integ.eri()
    nsh = basis.nshell()
    nao = basis.nbf()
    iνλσ = zeros(nao,nao,nao)
    Ci = C[:,i]
    for s in 1:nsh
        for r in 1:nsh
            for q in 1:nsh
                for p in 1:nsh
                    sf = basis.shell_to_basis_function(s-1) + 1
                    rf = basis.shell_to_basis_function(r-1) + 1
                    qf = basis.shell_to_basis_function(q-1) + 1
                    pf = basis.shell_to_basis_function(p-1) + 1
                    shell = mints.ao_eri_shell(p-1,q-1,r-1,s-1).np
                    pn, qn, rn, sn = size(shell)
                    pn -= 1
                    qn -= 1
                    rn -= 1
                    sn -= 1
                    #println("$pf $pn $qf $qn $rf $rn $sf $sn")
                    _Ci = Ci[pf:pf+pn]
                    @tensoropt begin
                        temp[ν,λ,σ] := _Ci[μ]*shell[μ,ν,λ,σ]
                    end
                    ##println("$p $q $r $s")
                    iνλσ[qf:qf+qn,rf:rf+rn,sf:sf+sn] += temp
                end
            end
        end
    end
    return iνλσ

    #ao_eri = DiskFourTensor("/tmp/disk_gao.$name.jues.0", Float64, nao, nao, nao, nao, "w")
    #blockfill!(ao_eri, 0.0)
    #while !si.is_done()
    #    p, q, r, s = (si.p, si.q, si.r, si.s)
    #    pf = basis.shell_to_basis_function(p) + 1
    #    qf = basis.shell_to_basis_function(q) + 1
    #    rf = basis.shell_to_basis_function(r) + 1
    #    sf = basis.shell_to_basis_function(s) + 1
    #    shell = mints.ao_eri_shell(p, q, r, s).to_array()
    #    pn, qn, rn, sn = size(shell)
    #    pn -= 1
    #    qn -= 1
    #    rn -= 1
    #    sn -= 1
    #    #ao_eri[pf:pf+pn, qf:qf+qn, sf:sf+sn, rf:rf+rn] = permutedims(shell, [1, 2, 4, 3])
    #    #ao_eri[qf:qf+qn, pf:pf+pn, sf:sf+sn, rf:rf+rn] =
    #    #    permutedims(permutedims(shell, [2, 1, 3, 4]), [1, 2, 4, 3])
    #    #ao_eri[rf:rf+rn, sf:sf+sn, pf:pf+pn, qf:qf+qn] =
    #    #    permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2])
    #    #ao_eri[sf:sf+sn, rf:rf+rn, pf:pf+pn, qf:qf+qn] = permutedims(
    #    #    permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2]),
    #    #    [2, 1, 3, 4],
    #    #)
    #    #ao_eri[rf:rf+rn, sf:sf+sn, qf:qf+qn, pf:pf+pn] = permutedims(
    #    #    permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2]),
    #    #    [1, 2, 4, 3],
    #    #)
    #    #ao_eri[sf:sf+sn, rf:rf+rn, qf:qf+qn, pf:pf+pn] =
    #    #    permutedims(permutedims(shell, [4, 2, 3, 1]), [1, 3, 2, 4])

    #    ao_eri1 = shell
    #    #ao_eri[pf:pf+pn, qf:qf+qn, rf:rf+rn, sf:sf+sn] = shell
    #    for ν in 1:qn
    #        νo = ν+qf
    #        for λ in 1:rn
    #            λo = λ+rf
    #            for σ in 1:sn
    #                σo = σ+sf
    #                tot = 0.0
    #                for μ in 1:pn
    #                    μo = μ+pf
    #                    tot += C[μo,i]*ao_eri1[μ,ν,λ,σ] 
    #                end
    #                iνλσ[νo,λo,σo] = tot
    #            end
    #        end
    #    end
    #    ao_eri2 = permutedims(shell, [2, 1, 3, 4])
    #    #ao_eri[qf:qf+qn, pf:pf+pn, rf:rf+rn, sf:sf+sn] = permutedims(shell, [2, 1, 3, 4])
    #    for \
    #    ao_eri3 = permutedims(shell, [1, 2, 4, 3])
    #    ao_eri4 = permutedims(permutedims(shell, [2, 1, 3, 4]), [1, 2, 4, 3])
    #    ao_eri5 = permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2])
    #    ao_eri6 = permutedims(
    #        permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2]),
    #        [2, 1, 3, 4],
    #    )
    #    ao_eri7 = permutedims(
    #        permutedims(permutedims(shell, [3, 2, 1, 4]), [1, 4, 3, 2]),
    #        [1, 2, 4, 3],
    #    )
    #    ao_eri8 =
    #        permutedims(permutedims(shell, [4, 2, 3, 1]), [1, 3, 2, 4])
    #    si.next()
    #end
end
