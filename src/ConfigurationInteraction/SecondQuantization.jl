module SecondQuantization

export Determinant
export αexcitation_level
export βexcitation_level
export excitation_level
export annihilate
export create
export exclusive
export phase

struct Determinant
    α::Int
    β::Int
end

function Determinant(α::String, β::String)
    
    αint = parse(Int, reverse(α); base=2) 
    βint = parse(Int, reverse(β); base=2) 

    Determinant(αint, βint)
end

function αexcitation_level(D1::Determinant, D2::Determinant)

    αexc = count(i->(i=='1'), bitstring(D1.α ⊻ D2.α))

    return αexc/2
end

function βexcitation_level(D1::Determinant, D2::Determinant)

    βexc = count(i->(i=='1'), bitstring(D1.β ⊻ D2.β))

    return βexc/2
end

function excitation_level(D1::Determinant, D2::Determinant)

    return αexcitation_level(D1, D2) + βexcitation_level(D1,D2)
end

function annihilate(D::Determinant, orb::Int, spin::Char)

    if spin == 'α'
        if D.α & (1 << (orb-1)) == 0
            error("Annihilation error. Orbital $orb is not occupied")
        end

        # Determine sign
        l = 0
        i = 1
        while i < (1 << (orb-1))
            l += D.α & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newα = D.α ⊻ (1 << (orb-1))

        return (-1)^l, Determinant(newα, D.β)

    elseif spin =='β' 
        if D.β & (1 << (orb-1)) == 0
            error("Annihilation error. Orbital $orb is not occupied")
        end

        # Determine sign
        l = count(i->(i=='1'), bitstring(D.α))
        i = 1
        while i < (1 << (orb-1))
            l += D.β & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newβ = D.β ⊻ (1 << (orb-1))

        return (-1)^l, Determinant(D.α, newβ)
    end
end

function create(D::Determinant, orb::Int, spin::Char)

    if spin == 'α'
        if D.α & (1 << (orb-1)) ≠ 0
            error("Creation error. Orbital $orb is occupied")
        end

        # Determine sign
        l = 0
        i = 1
        while i < (1 << (orb-1))
            l += D.α & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newα = D.α | (1 << (orb-1))

        return (-1)^l, Determinant(newα, D.β)

    elseif spin =='β' 
        if D.β & (1 << (orb-1)) ≠ 0
            error("Annihilation error. Orbital $orb is not occupied")
        end

        # Determine sign
        l = count(i->(i=='1'), bitstring(D.α))
        i = 1
        while i < (1 << (orb-1))
            l += D.β & i ≠ 0 ? 1 : 0
            i = i << 1
        end

        newβ = D.β | (1 << (orb-1))

        return (-1)^l, Determinant(D.α, newβ)
    end
end

function exclusive(D1::Determinant, D2::Determinant)

    αexcl = D1.α ⊻ D2.α & D1.α
    βexcl = D1.β ⊻ D2.β & D1.β

    out = []
    i = 1
    while 1<<(i-1) ≤ max(αexcl, βexcl)
        if 1<<(i-1) & αexcl ≠ 0
            push!(out, (i, 'α'))
        end
        if 1<<(i-1) & βexcl ≠ 0
            push!(out, (i, 'β'))
        end
        i += 1
    end
    return out
end

function phase(D1::Determinant, D2::Determinant)

    p = 1
    _det = Determinant(D1.α, D1.β)
    for (i,σ) in reverse(exclusive(D1, D2))
        f, _det = annihilate(_det, i, σ)
        p = f*p
    end
    for (a,σ) in exclusive(D1, D2)
        f, _det = create(_det, a, σ)
        p = f*p
    end

    return p
end

end
