module MollerPlesset
using JuES.Wavefunction
using JuES.DiskTensors
export do_rmp2
export do_ump2

function do_rmp2(refWfn::Wfn)
    dmp2 = 0.0
    nocc = refWfn.nalpha
    rocc = 1:1:refWfn.nalpha
    rvir = nocc+1:1:nocc+refWfn.nvira
#    @views moeri = permutedims(refWfn.pqrs,[1,3,2,4])
    moeri = permutedims(refWfn.ijab,[1,3,2,4])
    epsa = refWfn.epsa
    for b in rvir
    	for a in rvir
    	    for j in rocc
                for i in rocc
                    aa = a - nocc
                    bb = b - nocc
    	            dmp2 += (moeri[i,j,aa,bb]*(2*moeri[i,j,aa,bb] - moeri[i,j,bb,aa]))/
    		            (epsa[i] + epsa[j] - epsa[a] - epsa[b])
                end
    	    end
    	end
    end
    return dmp2
end
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
    #if typeof(refWfn.ijab) == Array{Float64,4}
    #    moeri = permutedims(refWfn.ijab,[1,3,2,4]) -
    #    permutedims(refWfn.ijab,[1,3,4,2])
    #end
    moeri = refWfn.ijab
    for b in rvira
    	for a in b+1:refWfn.nmo
            aa = a - nocca
            bb = b - nocca
            cache = moeri[:,aa:aa,:,bb:bb] - moeri[:,bb:bb,:,aa:aa]
            cache = permutedims(cache,[1,3,2,4])
    	    for j in rocca
    	        for i in j+1:refWfn.nalpha
                #dmp2 += (1/4)*(moeri[i,aa,j,bb] - moeri[i,bb,j,aa])*(moeri[i,aa,j,bb] - moeri[i,bb,j,aa])/
    	        #    (epsa[i] + epsa[j] - epsa[a] - epsa[b])
                moint = cache[i,j,1,1]
                dmp2 += (moint*moint)/
    	            (epsa[i] + epsa[j] - epsa[a] - epsa[b])
    	        end
    	    end
    	end
    end
    #spin case ABAB
    #if !unr
    #	moeri = permutedims(refWfn.ijab,[1,3,2,4])
    #else
    #	moeri = permutedims(refWfn.iJaB,[1,3,2,4])
    #end
    moeri = refWfn.iJaB
    for b in rvirb
    	for a in rvira
            aa = a - nocca
            bb = b - noccb
            cache = moeri[:,aa:aa,:,bb:bb]
            cache = permutedims(cache,[1,3,2,4])
    	    for j in roccb
    	        for i in rocca
                    dmp2 += cache[i,j,1,1]*cache[i,j,1,1]/(epsa[i] + epsb[j] - epsa[a] - epsb[b])
    	        end
    	    end
    	end
    end
    #spin case BBBB
    #if !unr
    #	moeri = permutedims(refWfn.ijab,[1,3,2,4]) - 
    #		permutedims(refWfn.ijab,[1,3,4,2])
    #else
    #	moeri = permutedims(refWfn.IJAB,[1,3,2,4]) - 
    #		permutedims(refWfn.IJAB,[1,3,4,2])
    #end
    moeri = refWfn.IJAB
    cache = zeros(noccb,1,noccb,1)
    for b in rvirb
    	for a in b+1:refWfn.nmo
            aa = a - noccb
            bb = b - noccb
            cache = moeri[:,aa:aa,:,bb:bb] .- moeri[:,bb:bb,:,aa:aa]
            cache = permutedims(cache,[1,3,2,4])
    	    for i in roccb
    	        for j in i+1:refWfn.nbeta
                    dmp2 += (cache[i,j,1,1]*cache[i,j,1,1])/
    		        (epsb[i] + epsb[j] - epsb[a] - epsb[b])
    		end
    	    end
    	end
    end
    return dmp2
    end
end #module
