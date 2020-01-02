module MollerPlesset
using Wavefunction
export do_rmp2
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
function transform_tei2(gao::Array{Float64,4},C::Array{Float64,2})
    norb = size(gao)[1]::Int64 #indexed from 1
    g1 = zeros(size(G))::Array{Float64,4}
    g2 = zeros(size(G))::Array{Float64,4}
    R = collect(UnitRange(1,norb))::Array{Int64,1}
    @inbounds @fastmath for s in R
        @inbounds @fastmath for si in R
            @inbounds @fastmath for rh in R
                @inbounds @fastmath for nu in R
                    @inbounds @fastmath @simd for mu in R
                        @views g1[mu,nu,rh,s] += gao[mu,nu,rh,si]*C[si,s]
                    end
                end
            end
        end
    end
    @inbounds @fastmath for s in R
        @inbounds @fastmath for r in R
            @inbounds @fastmath for rh in R
                @inbounds @fastmath for nu in R
                    @inbounds @fastmath @simd for mu in R
                        @views g2[mu,nu,r,s] += g1[mu,nu,rh,s]*C[rh,r]
                    end
                end
            end
        end
    end
    g1 .= 0.0
    @inbounds @fastmath for s in R
        @inbounds @fastmath for r in R
            @inbounds @fastmath for q in R
                @inbounds @fastmath for nu in R
                    @inbounds @fastmath @simd for mu in R
                        @views g1[mu,q,r,s] += g2[mu,nu,r,s]*C[nu,q]
                    end
                end
            end
        end
    end
    g2 .= 0.0
    @inbounds @fastmath for s in R
        @inbounds @fastmath for r in R
            @inbounds @fastmath for q in R
                @inbounds @fastmath for mu in R
                    @inbounds @fastmath @simd for p in R
                        @views g2[p,q,r,s] += g1[mu,q,r,s]*C[mu,p]
                    end
                end
            end
        end
    end
    return g2
end
#function mp2(moeri::Array{Float64,4},e::Array{Float64,1},nocc::Int64)
#    moeri_alt = permutedims(moeri,[1,2,4,3])::Array{Float64,4}
#    delta_mp2 = 0.0::Float64
#    for I in 1:1:nocc
#        for J in 1:1:nocc
#            for A in nocc+1:1:length(e)
#                for B in nocc+1:1:length(e)
#                    delta_mp2 += (moeri[I,J,A,B]*(2*moeri[I,J,A,B] - moeri_alt[I,J,A,B]))/(e[I] + e[J] - e[A] - e[B])
#                end
#            end
#        end
#    end
#    return delta_mp2
#end

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
