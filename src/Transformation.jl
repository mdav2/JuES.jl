module Transformation
using JuES.DiskTensors
export disk_tei_transform
export mem_tei_transform

function disk_tei_transform(gao::Array{Float64,4},
                       C1::Array{Float64,2},
                       C2::Array{Float64,2},
                       C3::Array{Float64,2},
                       C4::Array{Float64,2},
                       name::String)
    """
    disk based, general spin orbital transformation
    """
    norb = size(C1)[1]
    rr1 = UnitRange(1,size(C1)[2])
    rr2 = UnitRange(1,size(C2)[2])
    rr3 = UnitRange(1,size(C3)[2])
    rr4 = UnitRange(1,size(C4)[2])
    R1 = collect(rr1)
    R2 = collect(rr2)
    R3 = collect(rr3)
    R4 = collect(rr4)
    g1 = DiskFourTensor("/tmp/jues.$name.g1.0",Float64,norb,norb,norb,norb,"w")
    g2 = DiskFourTensor("/tmp/jues.$name.0",Float64,norb,norb,norb,norb,"w")
    blockfill!(g1,0.0)
    blockfill!(g2,0.0)
    cache = zeros(norb,norb)
    cache2 = zeros(norb,norb)
    for s in R
    	for rh in R
    	    for si in R
    	        cache[:,:] = gao[rr,rr,rh,si]
    		    for nu in R
                    for mu in R
                        cache2[mu,nu] += cache[mu,nu]*C1[si,s]
                    end
                end
            end
            g1[rr,rr,rh,s] = cache2[:,:]
            cache2[:,:] = zeros(norb,norb)
        end
    end
    for s in R
        for r in R
            for rh in R
                cache[:,:] = g1[:,:,rh,s]
                for nu in R
                    for mu in R
                        #@views g2[mu,nu,r,s] += g1[mu,nu,rh,s]*C[rh,r]
			@views cache2[mu,nu] += cache[mu,nu]*C2[rh,r]
                    end
                end
            end
            g2[:,:,r,s] = cache2
            cache2[:,:] = zeros(norb,norb)
        end
    end
    g1[:,:,:,:] = 0.0
    for s in R
        for r in R
	    cache = g2[:,:,r,s]
            #cache2[:,:] = 0.0
            for q in R
                for nu in R
                    for mu in R
                        #@views g1[mu,q,r,s] += g2[mu,nu,r,s]*C[nu,q]
                        cache2[mu,q] += cache[mu,nu]*C3[nu,q]
                    end
                end
            end
            g1[:,:,r,s] = cache2
            cache2[:,:] = zeros(norb,norb)
        end
    end
    g2[:,:,:,:] = 0.0
    for s in R
        for r in R
            cache = g1[:,:,r,s]
            for q in R
                for p in R
                    for mu in R
                        #@views g2[p,q,r,s] += g1[mu,q,r,s]*C[mu,p]
			cache2[p,q] += cache[mu,q]*C4[mu,p]
                    end
                end
            end
            g2[:,:,r,s] = cache2
            cache2[:,:] = zeros(norb,norb)
        end
    end
    return g2
end


"""
disk based, restricted transformation
"""
function disk_tei_transform(gao::Array{Float64,4},C::Array{Float64,2},name::String)
    norb = size(C)[1]
    rr = UnitRange(1,norb)
    R = collect(UnitRange(1,norb))::Array{Int64,1}
    g1 = DiskFourTensor("/tmp/jues.$name.g1.0",Float64,norb,norb,norb,norb,"w")
    g2 = DiskFourTensor("/tmp/jues.$name.0",Float64,norb,norb,norb,norb,"w")
    blockfill!(g1,0.0)
    blockfill!(g2,0.0)
    cache = zeros(norb,norb)
    cache2 = zeros(norb,norb)
    for s in R
        for rh in R
            for si in R
                cache[:,:] = gao[rr,rr,rh,si]
                for nu in R
                    for mu in R
                        #g1[mu,nu,rh,s] += gao[mu,nu,rh,si]*C[si,s]
                        cache2[mu,nu] += cache[mu,nu]*C[si,s]
                    end
                end
            end
	    g1[rr,rr,rh,s] = cache2[:,:]
	    cache2[:,:] = zeros(norb,norb)
        end
    end
    for s in R
        for r in R
            for rh in R
		        cache[:,:] = g1[:,:,rh,s]
                for nu in R
                    for mu in R
                        #@views g2[mu,nu,r,s] += g1[mu,nu,rh,s]*C[rh,r]
			            cache2[mu,nu] += cache[mu,nu]*C[rh,r]
                    end
                end
            end
	    g2[:,:,r,s] = cache2
	    cache2[:,:] = zeros(norb,norb)
        end
    end
    g1[:,:,:,:] = 0.0
    for s in R
        for r in R
	    cache = g2[:,:,r,s]
	    #cache2[:,:] = 0.0
            for q in R
                for nu in R
                    for mu in R
                        #@views g1[mu,q,r,s] += g2[mu,nu,r,s]*C[nu,q]
			            cache2[mu,q] += cache[mu,nu]*C[nu,q]
                    end
                end
            end
	    g1[:,:,r,s] = cache2
	    cache2[:,:] = zeros(norb,norb)
        end
    end
    g2[:,:,:,:] = 0.0
    for s in R
        for r in R
	    cache = g1[:,:,r,s]
            for q in R
                for p in R
                    for mu in R
                        #@views g2[p,q,r,s] += g1[mu,q,r,s]*C[mu,p]
			            cache2[p,q] += cache[mu,q]*C[mu,p]
                    end
                end
            end
	        g2[:,:,r,s] = cache2
            cache2[:,:] = zeros(norb,norb)
        end
    end
    return g2
end
function mem_tei_transform(gao::Array{Float64,4},C::Array{Float64,2})
    norb = size(gao)[1]::Int64 #indexed from 1
    g1 = zeros(size(gao))::Array{Float64,4}
    g2 = zeros(size(gao))::Array{Float64,4}
    R = collect(UnitRange(1,norb))::Array{Int64,1}
    for s in R
        for si in R
            for rh in R
                for nu in R
                    @simd for mu in R
                        @views g1[mu,nu,rh,s] += gao[mu,nu,rh,si]*C[si,s]
                    end
                end
            end
        end
    end
    for s in R
        for r in R
            for rh in R
                for nu in R
                    @simd for mu in R
                        @views g2[mu,nu,r,s] += g1[mu,nu,rh,s]*C[rh,r]
                    end
                end
            end
        end
    end
    g1 .= 0.0
    for s in R
        for r in R
            for q in R
                for nu in R
                    @simd for mu in R
                        @views g1[mu,q,r,s] += g2[mu,nu,r,s]*C[nu,q]
                    end
                end
            end
        end
    end
    g2 .= 0.0
    for s in R
        for r in R
            for q in R
                for mu in R
                    @simd for p in R
                        @views g2[p,q,r,s] += g1[mu,q,r,s]*C[mu,p]
                    end
                end
            end
        end
    end
    return g2
end
end
