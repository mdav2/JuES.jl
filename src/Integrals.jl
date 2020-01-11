module Integrals

using JuES.DiskTensors
using PyCall
const psi4 = PyNULL()
function __init__()
    copy!(psi4, pyimport("psi4"))
end

export disk_ao

"""
    disk_ao

computes all AO basis TEI and fills a DiskFourTensor object with those
"""
function disk_ao(mints::PyObject, basis::PyObject, name::String = "default")
    integ = mints.integral()
    si = integ.shells_iterator()
    si.first()
    nao = basis.nao()
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

end #module
