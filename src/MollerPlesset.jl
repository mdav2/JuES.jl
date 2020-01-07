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
    moeri = permutedims(refWfn.ijab,[1,3,2,4]) -
    permutedims(refWfn.ijab,[1,3,4,2])
    for b in rvira
    	for a in rvira
    	    for j in rocca
    		    for i in rocca
                    aa = a - nocca
                    bb = b - nocca
    		        dmp2 += (1/4)*moeri[i,j,aa,bb]*moeri[i,j,aa,bb]/
    		        (epsa[i] + epsa[j] - epsa[a] - epsa[b])
    		    end
    	    end
    	end
    end
    #spin case ABAB
    if !unr
    	moeri = permutedims(refWfn.ijab,[1,3,2,4])
    else
    	moeri = permutedims(refWfn.iJaB,[1,3,2,4])
    end
    for b in rvirb
    	for a in rvira
    	    for j in roccb
    		    for i in rocca
                    aa = a - nocca
                    bb = b - noccb
                    dmp2 += moeri[i,j,aa,bb]*(moeri[i,j,aa,bb]/(epsa[i] + epsb[j] - epsa[a] - epsb[b]))
    		    end
    	    end
    	end
    end
    #spin case BBBB
    if !unr
    	moeri = permutedims(refWfn.ijab,[1,3,2,4]) - 
    		permutedims(refWfn.ijab,[1,3,4,2])
    else
    	moeri = permutedims(refWfn.IJAB,[1,3,2,4]) - 
    		permutedims(refWfn.IJAB,[1,3,4,2])
    end
    for b in rvirb
    	for a in rvirb
    	    for i in roccb
    	        for j in roccb
                    aa = a - noccb
                    bb = b - noccb
                    dmp2 += (1/4)*moeri[i,j,aa,bb]*moeri[i,j,aa,bb]/
    		        (epsb[i] + epsb[j] - epsb[a] - epsb[b])
    		    end
    	    end
    	end
    end
    return dmp2
    end
end #module
