"""
    do_ump2

performs unrestricted Moller-Plesset perturbation theory to second order
on a UHF reference.
## paramters
refWfn::Wfn -> wavefunction to which the MP2 equations will be applied.

## outputs
dmp2::Float -> UMP2 energy correction
"""
function do_ump2(refWfn::Wfn)
    dmp2 = 0.0
    unr = refWfn.unrestricted
    nocca = refWfn.nalpha
    noccb = refWfn.nbeta
    rocca = 1:refWfn.nalpha
    roccb = 1:refWfn.nbeta
    rvira = refWfn.nalpha+1:refWfn.nmo
    rvirb = refWfn.nbeta+1:refWfn.nmo
    epsa = refWfn.epsa
    epsb = refWfn.epsb
    #spin case AAAA
    moeri = refWfn.ijab
    for b in rvira
        for a = b+1:refWfn.nmo
            aa = a - nocca
            bb = b - nocca
            cache = moeri[:, aa:aa, :, bb:bb] - moeri[:, bb:bb, :, aa:aa]
            cache = squeeze(cache)
            for j in rocca
                for i = j+1:refWfn.nalpha
                    moint = cache[i, j]
                    dmp2 += (moint * moint) / (epsa[i] + epsa[j] - epsa[a] - epsa[b])
                end
            end
        end
    end
    #spin case ABAB
    moeri = refWfn.iJaB
    for b in rvirb
        for a in rvira
            aa = a - nocca
            bb = b - noccb
            cache = moeri[:, aa:aa, :, bb:bb]
            cache = squeeze(cache)
            for j in roccb
                for i in rocca
                    dmp2 +=
                        cache[i, j] * cache[i, j] / (epsa[i] + epsb[j] - epsa[a] - epsb[b])
                end
            end
        end
    end
    #spin case BBBB
    moeri = refWfn.IJAB
    cache = zeros(noccb, 1, noccb, 1)
    for b in rvirb
        for a = b+1:refWfn.nmo
            aa = a - noccb
            bb = b - noccb
            cache = moeri[:, aa:aa, :, bb:bb] .- moeri[:, bb:bb, :, aa:aa]
            cache = squeeze(cache)
            for i in roccb
                for j = i+1:refWfn.nbeta
                    dmp2 +=
                        (cache[i, j] * cache[i, j]) /
                        (epsb[i] + epsb[j] - epsb[a] - epsb[b])
                end
            end
        end
    end
    return dmp2
end
