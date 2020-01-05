module MollerPlesset
using JuES.Wavefunction
using JuES.DiskTensors
export do_rmp2
export transform_tei
export transform_tei2
function gaogen(refWfn::Wfn)

end
function direct_MO(gaogen,Wfn,p,q,r,s)
	R = 1:1:nmo
	g2 = zeros(nmo,nmo)
	for si in R
		for rh in R
			for nu in R
				for mu in R
					tsum = 0.0
					for si in R
						tsum += gaogen(mu,nu,rh,si)*C[si,s]
					end
					g2[mu,nu] += tsum*C[rh,r]
				end
			end
		end
	end
end
function transform_tei(gao::Array{Float64,4},C::Array{Float64,2})
	# doing transform by matrix 
	# specify p and q
	# so all r and all s
	norb = size(gao)[1]
	rr = UnitRange(1,norb)
    R = collect(UnitRange(1,norb))::Array{Int64,1}
	g1 = DiskFourTensor("/tmp/g1.h5",Float64,norb,norb,norb,norb,"w")
	g2 = DiskFourTensor("/tmp/g2.h5",Float64,norb,norb,norb,norb,"w")
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
						@views cache2[mu,nu] += cache[mu,nu]*C[rh,r]
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
function transform_tei2(gao::Array{Float64,4},C::Array{Float64,2})
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

function do_rmp2(refWfn::Wfn)
	dmp2 = 0.0
	nocc = refWfn.nalpha
	rocc = 1:1:refWfn.nalpha
	rvir = nocc+1:1:nocc+refWfn.nvira
	@views moeri = permutedims(refWfn.pqrs,[1,3,2,4])
	epsa = refWfn.epsa
	for b in rvir
		for a in rvir
			for j in rocc
				for i in rocc
					dmp2 += (moeri[i,j,a,b]*(2*moeri[i,j,a,b] - moeri[i,j,b,a]))/
							 (epsa[i] + epsa[j] - epsa[a] - epsa[b])
				end
			end
		end
	end
	return dmp2
end
end #module
