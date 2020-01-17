module Transformation
using JuES.DiskTensors
using TensorOperations
using Base.Threads
using LinearAlgebra
export disk_tei_transform
export mem_tei_transform
export tei_transform

function tei_transform(
    gao::Array{Float64,4},
    C1::Array{Float64,2},
    C2::Array{Float64,2},
    C3::Array{Float64,2},
    C4::Array{Float64,2},
    name::String,
)
    """
    general spin orbital transformation
     
    """
    T = typeof(gao)
    norb = size(C1)[1]
    d = size(C1)[1]
    d1 = size(C1)[2]
    d2 = size(C2)[2]
    d3 = size(C3)[2]
    d4 = size(C4)[2]
    rr = UnitRange(1, norb)
    rr1 = UnitRange(1, d1)
    rr2 = UnitRange(1, d2)
    rr3 = UnitRange(1, d3)
    rr4 = UnitRange(1, d4)
    R = collect(rr)
    R1 = collect(rr1)
    R2 = collect(rr2)
    R3 = collect(rr3)
    R4 = collect(rr4)
    temp = zeros(d, d, d, d4)
    #Quarter transform 1
    #(μ,ν,λ,σ) -> (μ,ν,λ,b)
    @tensor begin
        temp[μ,ν,λ,b] = C4[σ,b]*gao[μ,ν,λ,σ]
    end
    #for b in R4
    #    for μ in R
    #        ocache = zeros(d, d)
    #        for ν in R
    #            icache = gao[μ, ν, :, :]
    #            for λ in R
    #                #@views ocache[ν,λ] = dot(C4[:,b],icache[λ,:])
    #                for σ in R
    #                    ocache[ν, λ] += C4[σ, b] * icache[λ, σ]
    #                end
    #            end
    #        end
    #        temp[μ, :, :, b] = ocache[:, :]
    #    end
    #end
    #Quarter transform 2
    #(μ,ν,λ,b) -> (μ,ν,j,b)
    temp2 = zeros(d, d, d3, d4)
    @tensor begin
        temp2[μ,ν,j,b] = C3[λ,j]*temp[μ,ν,λ,b]
    end
    #for b in R4
    #    for j in R3
    #        ocache = zeros(d, d)
    #        for μ in R
    #            icache = temp[μ, :, :, b]
    #            for ν in R
    #                for λ in R #contract over λ
    #                    ocache[μ, ν] += C3[λ, j] * icache[ν, λ]
    #                end
    #            end
    #        end
    #        temp2[:, :, j, b] = ocache[:, :]
    #    end
    #end
    #Quarter transform 3
    #(μ,ν,j,b) -> (μ,a,j,b)
    temp = zeros(d, d2, d3, d4)
    @tensor begin
        temp[μ,a,j,b] = C2[ν,a]*temp2[μ,ν,j,b]
    end
    #for b in R4
    #    for j in R3
    #        ocache = zeros(d, d2)
    #        icache = temp2[:, :, j, b]
    #        for a in R2
    #            for μ in R
    #                for ν in R
    #                    #temp[μ,a,j,b] += C2[ν,a]*temp2[μ,ν,j,b]
    #                    ocache[μ, a] += C2[ν, a] * icache[μ, ν]
    #                end
    #            end
    #        end
    #        temp[:, :, j, b] = ocache[:, :]
    #    end
    #end
    #Quarter transform 4
    #(μ,a,j,b) -> (i,a,j,b)

    temp2 = zeros(d1, d2, d3, d4)
    @tensor begin
        temp2[i,a,j,b] = C1[μ,i]*temp[μ,a,j,b]
    end
    #for b in R4
    #    for j in R3
    #        icache = temp[:, :, j, b]
    #        ocache = zeros(d1, d2)
    #        for a in R2
    #            for i in R1
    #                for μ in R
    #                    #temp2[i,a,j,b] += C1[μ,i]*temp[μ,a,j,b]
    #                    ocache[i, a] += C1[μ, i] * icache[μ, a]
    #                end
    #            end
    #        end
    #        temp2[:, :, j, b] = ocache
    #    end
    #end
    return temp2
end
function tei_transform(
    gao::DiskFourTensor,
    C1::Array{Float64,2},
    C2::Array{Float64,2},
    C3::Array{Float64,2},
    C4::Array{Float64,2},
    name::String,
)
    """
    general spin orbital transformation
     
    """
    T = typeof(gao)
    norb = size(C1)[1]
    d = size(C1)[1]
    d1 = size(C1)[2]
    d2 = size(C2)[2]
    d3 = size(C3)[2]
    d4 = size(C4)[2]
    rr = UnitRange(1, norb)
    rr1 = UnitRange(1, d1)
    rr2 = UnitRange(1, d2)
    rr3 = UnitRange(1, d3)
    rr4 = UnitRange(1, d4)
    R = collect(rr)
    R1 = collect(rr1)
    R2 = collect(rr2)
    R3 = collect(rr3)
    R4 = collect(rr4)
    if T == Array{Float64,4}
        temp = zeros(d, d, d, d4)
    elseif T == DiskFourTensor
        temp = DiskFourTensor("/tmp/jues.$name.temp.0", Float64, d, d, d, d4, "w")
        blockfill!(temp, 0.0)
    end
    #Quarter transform 1
    #(μ,ν,λ,σ) -> (μ,ν,λ,b)
    #for μ in R
    #    for ν in R
    #        ocache = zeros(d,d4)
    #        icache = gao[μ,ν,:,:]
    #        @tensoropt begin
    #            ocache[λ,b] = C4[σ,b]*icache[λ,σ]
    #        end
    #        temp[μ,ν,:,:] = ocache[:,:]
    #    end
    #end
    for μ in R
        for ν in R
            ocache = zeros(d,d4)
            icache = gao[μ,ν,:,:]
            for λ in R 
                for b in R4
                    for σ in R
                        ocache[λ,b] += C4[σ,b]*icache[λ,σ]
                    end
                end
            end
            #@tensoropt begin
            #    ocache[λ,b] = C4[σ,b]*icache[λ,σ]
            #end
            temp[μ,ν,:,:] = ocache[:,:]
        end
    end
    #for b in R4
    #    for μ in R
    #        ocache = zeros(d, d)
    #        for ν in R
    #            icache = gao[μ, ν, :, :]
    #            temp = zeros(d)
    #            @tensoropt begin
    #                temp[λ] = C4[σ,b]*icache[λ,σ]
    #            end
    #            ocache[ν,:] = temp[:]
    #            #for λ in R
    #            #    for σ in R
    #            #        ocache[ν, λ] += C4[σ, b] * icache[λ, σ]
    #            #    end
    #            #end
    #        end
    #        temp[μ, :, :, b] = ocache[:, :]
    #    end
    #end
    #Quarter transform 2
    #(μ,ν,λ,b) -> (μ,ν,j,b)
    if T == Array{Float64,4}
        temp2 = zeros(d, d, d3, d4)
    elseif T == DiskFourTensor
        temp2 = DiskFourTensor("/tmp/jues.$name.temp2.0", Float64, d, d, d3, d4, "w")
        blockfill!(temp2, 0.0)
    end
    for b in R4
        for j in R3
            ocache = zeros(d, d)
            for μ in R
                icache = temp[μ, :, :, b]
                for ν in R
                    for λ in R #contract over λ
                        ocache[μ, ν] += C3[λ, j] * icache[ν, λ]
                    end
                end
            end
            temp2[:, :, j, b] = ocache[:, :]
        end
    end
    #Quarter transform 3
    #(μ,ν,j,b) -> (μ,a,j,b)
    if T == Array{Float64,4}
        temp = zeros(d, d2, d3, d4)
    elseif T == DiskFourTensor
        temp = DiskFourTensor("/tmp/jues.$name.temp.0", Float64, d, d2, d3, d4, "w")
        blockfill!(temp, 0.0)
    end
    for b in R4
        for j in R3
            ocache = zeros(d, d2)
            icache = temp2[:, :, j, b]
            for a in R2
                for μ in R
                    for ν in R
                        #temp[μ,a,j,b] += C2[ν,a]*temp2[μ,ν,j,b]
                        ocache[μ, a] += C2[ν, a] * icache[μ, ν]
                    end
                end
            end
            temp[:, :, j, b] = ocache[:, :]
        end
    end
    #Quarter transform 4
    #(μ,a,j,b) -> (i,a,j,b)

    if T == Array{Float64,4}
        temp2 = zeros(d1, d2, d3, d4)
    elseif T == DiskFourTensor
        temp2 = DiskFourTensor("/tmp/jues.$name.temp2.0", Float64, d1, d2, d3, d4, "w")
        blockfill!(temp2, 0.0)
    end
    for b in R4
        for j in R3
            icache = temp[:, :, j, b]
            ocache = zeros(d1, d2)
            for a in R2
                for i in R1
                    for μ in R
                        #temp2[i,a,j,b] += C1[μ,i]*temp[μ,a,j,b]
                        ocache[i, a] += C1[μ, i] * icache[μ, a]
                    end
                end
            end
            temp2[:, :, j, b] = ocache
        end
    end
    return temp2
end


"""
disk based, restricted transformation
"""
function disk_tei_transform(gao::Array{Float64,4}, C::Array{Float64,2}, name::String)
    norb = size(C)[1]
    rr = UnitRange(1, norb)
    R = collect(UnitRange(1, norb))::Array{Int64,1}
    g1 = DiskFourTensor("/tmp/jues.$name.g1.0", Float64, norb, norb, norb, norb, "w")
    g2 = DiskFourTensor("/tmp/jues.$name.0", Float64, norb, norb, norb, norb, "w")
    blockfill!(g1, 0.0)
    blockfill!(g2, 0.0)
    cache = zeros(norb, norb)
    cache2 = zeros(norb, norb)
    for s in R
        for rh in R
            for si in R
                cache[:, :] = gao[rr, rr, rh, si]
                for nu in R
                    for mu in R
                        #g1[mu,nu,rh,s] += gao[mu,nu,rh,si]*C[si,s]
                        cache2[mu, nu] += cache[mu, nu] * C[si, s]
                    end
                end
            end
            g1[rr, rr, rh, s] = cache2[:, :]
            cache2[:, :] = zeros(norb, norb)
        end
    end
    for s in R
        for r in R
            for rh in R
                cache[:, :] = g1[:, :, rh, s]
                for nu in R
                    for mu in R
                        #@views g2[mu,nu,r,s] += g1[mu,nu,rh,s]*C[rh,r]
                        cache2[mu, nu] += cache[mu, nu] * C[rh, r]
                    end
                end
            end
            g2[:, :, r, s] = cache2
            cache2[:, :] = zeros(norb, norb)
        end
    end
    g1[:, :, :, :] = 0.0
    for s in R
        for r in R
            cache = g2[:, :, r, s]
            #cache2[:,:] = 0.0
            for q in R
                for nu in R
                    for mu in R
                        #@views g1[mu,q,r,s] += g2[mu,nu,r,s]*C[nu,q]
                        cache2[mu, q] += cache[mu, nu] * C[nu, q]
                    end
                end
            end
            g1[:, :, r, s] = cache2
            cache2[:, :] = zeros(norb, norb)
        end
    end
    g2[:, :, :, :] = 0.0
    for s in R
        for r in R
            cache = g1[:, :, r, s]
            for q in R
                for p in R
                    for mu in R
                        #@views g2[p,q,r,s] += g1[mu,q,r,s]*C[mu,p]
                        cache2[p, q] += cache[mu, q] * C[mu, p]
                    end
                end
            end
            g2[:, :, r, s] = cache2
            cache2[:, :] = zeros(norb, norb)
        end
    end
    return g2
end
function mem_tei_transform(gao::Array{Float64,4}, C::Array{Float64,2})
    norb = size(gao)[1]::Int64 #indexed from 1
    g1 = zeros(size(gao))::Array{Float64,4}
    g2 = zeros(size(gao))::Array{Float64,4}
    R = collect(UnitRange(1, norb))::Array{Int64,1}
    for s in R
        for si in R
            for rh in R
                for nu in R
                    @simd for mu in R
                        @views g1[mu, nu, rh, s] += gao[mu, nu, rh, si] * C[si, s]
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
                        @views g2[mu, nu, r, s] += g1[mu, nu, rh, s] * C[rh, r]
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
                        @views g1[mu, q, r, s] += g2[mu, nu, r, s] * C[nu, q]
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
                        @views g2[p, q, r, s] += g1[mu, q, r, s] * C[mu, p]
                    end
                end
            end
        end
    end
    return g2
end
end
