"""
    do_rmp2

performs MP2 on a restricted HF reference.
## paramters
refWfn::Wfn -> wavefunction to which MP2 will be applied.

## outputs
dmp2::Float MP2 energy correction
"""
function do_rmp2(refWfn::Wfn)
    dmp2 = 0.0
    nocc = refWfn.nalpha
    rocc = 1:1:refWfn.nalpha
    rvir = nocc+1:1:nocc+refWfn.nvira
    #    @views moeri = permutedims(refWfn.pqrs,[1,3,2,4])
    #
    Cao = refWfn.Cao
    Cav = refWfn.Cav
    moeri = permutedims(tei_transform(refWfn.uvsr, Cao, Cav, Cao, Cav, "oovv"), [1, 3, 2, 4])
    epsa = refWfn.epsa
    for b in rvir
        for a in rvir
            for j in rocc
                for i in rocc
                    aa = a - nocc
                    bb = b - nocc
                    dmp2 +=
                        (
                            moeri[i, j, aa, bb] *
                            (2 * moeri[i, j, aa, bb] - moeri[i, j, bb, aa])
                        ) / (epsa[i] + epsa[j] - epsa[a] - epsa[b])
                end
            end
        end
    end
    return dmp2
end
