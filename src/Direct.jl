module Direct
using PyCall
using JuES.Wavefunction
const psi4 = PyNULL()
function __init__()
    copy!(psi4, pyimport("psi4"))
end
export ao
export func_shell
export ao_to_mo_shell

function ao_function(wfn::DirectWfn, p, q, r, s)
    """
    Computes a single ERI in AO basis
    Very wasteful (throws out rest of shell)
    use for debugging
    """
    shell = ao_shell(wfn, p, q, r, s)
    return shell[po, qo, ro, so]
end
function ao_to_mo_shell(p, q, r, s, C, wfn::DirectWfn)
    nao = size(C)[1]
    moint = 0.0
    nshell = wfn.basis.nshell()
    aoshell = func_shell(wfn)
    println(aoshell)
    for shell1 = 1:nshell
        nao_shell1 = aoshell[shell1]
        for shell2 = 1:nshell
            nao_shell2 = aoshell[shell2]
            for shell3 = 1:nshell
                nao_shell3 = aoshell[shell3]
                for shell4 = 1:nshell
                    nao_shell4 = aoshell[shell4]
                    current = ao_shell(wfn, shell1, shell2, shell3, shell4)
                    for p in nao_shell1
                        for q in nao_shell2
                            for r in nao_shell3
                                for s in nao_shell4
                                    println(p, q, r, s)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
function ao_to_mo(p, q, r, s, C, wfn::DirectWfn)
    nao = size(C)[1]
    moint = 0.0
    for uu = 1:nao
        for vv = 1:nao
            for ss = 1:nao
                for rr = 1:nao
                    moint +=
                        ao_function(wfn, uu, vv, ss, rr) *
                        C[p, uu] *
                        C[q, vv] *
                        C[r, ss] *
                        C[s, rr]
                end
            end
        end
    end
    return moint
end
function func_shell(wfn::DirectWfn)
    """
    Matches basis function number and corresponding shell number
    usage: output[shellno] -> [fn1,fn2,...,fnn]
    """
    basis = wfn.basis
    nao = basis.nao()
    output = Dict()
    for i = 1:1:nao
        x = basis.function_to_shell(i - 1) + 1
        if haskey(output, x)
            push!(output[x], i)
        else
            output[x] = [i]
        end
    end
    println(output)
    return output
end
function ao_shell(wfn::DirectWfn, p, q, r, s)
    """
    Computes the AO eri shell (pq|rs)
    """
    basis = wfn.basis
    mints = wfn.mints
    ps = basis.function_to_shell(p - 1)
    qs = basis.function_to_shell(q - 1)
    rs = basis.function_to_shell(r - 1)
    ss = basis.function_to_shell(s - 1)

    po = p - basis.shell_to_ao_function(ps)# + 1
    qo = q - basis.shell_to_ao_function(qs)# + 1
    ro = r - basis.shell_to_ao_function(rs)# + 1
    so = s - basis.shell_to_ao_function(ss)# + 1
    shell = mints.ao_eri_shell(ps, qs, rs, ss).to_array()
    return shell
end
function direct_MP2(wfn::DirectWfn)
    """
    Algorithm
    --
    for each occ <ij|-->, compute partially transformed
    integrals <ij|λσ>, store in memory

    for each vir <--|ab> contract <ij|λσ> ⋅ C[a,λ] ⋅ C[b,σ]
    and its transpose <--|ba>

    compute contribution of <ij|ab> to δE(MP2)

       12 34      12 34     12 34        
    = <ij|ab> ( 2<ij|ab> - <ij|ba> ) / Δ[ijab]
    = <ij|ab> ( 2<ij|ab> - <ij|ba> ) / Δ[ijab]
    ------------------------------------------
       13 24      13 24     13 24        
    = [ia|jb] ( 2[ia|jb] - [ib|ja] ) / Δ[ijab]
    """
    basis = wfn.basis
    nshell = basis.nshell()
    #figure out how to make use of shell/batches
    # "         how to make partially transformed integrals
    # "         how to do contraction
    #for storing partially transformed integrals
    partial = zeros(Float64, wfn.nalpha, wfn.nalpha, wfn.nmo, wfn.nmo)
    R = collect(UnitRange(1, nshell))
    aoshell = func_shell(wfn)

end
end#module
