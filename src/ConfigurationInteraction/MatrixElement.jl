module MatrixElement
using TensorOperations
using JuES.ConfigurationInteraction.DetOperations

export Hd0
export Hd1
export Hd2

"""
Module implementing matrix elements (currently just Hamiltonian) via
Slater's rules

Formulas and notation from Szabo and Ostlund p. 70

syntax --> <nparticle>Electron_<matrix>d<ndiff>
	   -->
	   --> e.g. for the one particle matrix element of a Hamiltonian betweeen
	   --> two determinants differing by 1 orbital:
	   --> 	OneElectron_Hd1
"""

using JuES.Wavefunction
using JuES.ConfigurationInteraction.DetOperations

export Hd0
export Hd1
export Hd2

function Hd0(αindex::Array{Int64,1}, βindex::Array{Int64,1}, h::Array{Float64, 2}, V::Array{Float64, 4})
    """
    Σ [m|h|m] + 1/2 ΣΣ [mm|nn] - [mn|nm] 
    """

    # One-electron contribution
    _hαα = h[αindex, αindex]
    _hββ = h[βindex, βindex]
    @tensor E = _hαα[m,m] + _hββ[m,m]

    # Two-electron contributions

    _Vα  = V[αindex, αindex, αindex, αindex]
    _Vβ  = V[βindex, βindex, βindex, βindex]
    _Vαβ = V[αindex, αindex, βindex, βindex]

    @tensor E += 0.5*(_Vα[m,m,n,n] + _Vβ[m,m,n,n] + 2_Vαβ[m,m,n,n] - _Vα[m,n,n,m] - _Vβ[m,n,n,m])

    return E

end

function Hd1(αindex::Array{Int64,1}, βindex::Array{Int64,1}, D1::Determinant, D2::Determinant, h::Array{Float64,2}, V::Array{Float64, 4})
    """
    differ m -> p
    [m|h|p] + Σ[mp|nn] - [mn|np]
    """

    p = phase(D1, D2)

    # if m and p are α
    if αexcitation_level(D1, D2) == 1
        m, = αexclusive_index(D1, D2)
        p, = αexclusive_index(D2, D1)

        _Jα = V[m, p, αindex, αindex]
        _Jβ = V[m, p, βindex, βindex]
        _Kα = V[m, αindex, αindex, p]

        @tensor E = _Jα[n,n] + _Jβ[n,n] - _Kα[n,n]

        return p*(h[m,p] + E)

    else
        m, = βexclusive_index(D1, D2)
        p, = βexclusive_index(D2, D1)

        _Jα = V[m, p, αindex, αindex]
        _Jβ = V[m, p, βindex, βindex]
        _Kβ = V[m, βindex, βindex, p]

        @tensor E = _Jα[n,n] + _Jβ[n,n] - _Kβ[n,n]

        return p*(h[m,p] + E)
    end
end

function Hd2(D1::Determinant, D2::Determinant, V::Array{Float64, 4})
    """
    mn -> pq
    [mp|nq] - [mq|np]
    """

    p = phase(D1, D2)

    # If α excitation is one, it means m and n have different spins 
    if αexcitation_level(D1, D2) == 1
        m, = αexclusive_index(D1, D2)
        n, = βexclusive_index(D1, D2)
        p, = αexclusive_index(D2, D1)
        q, = βexclusive_index(D2, D1)

        if (D1.α, D1.β) == (1,1) || (D2.α, D2.β) == (1,1)
            println("\nHd2")
            showdet(D1, 2)
            showdet(D2, 2)
            println("---")
            println(m)
            println(n)
            println(p)
            println(q)
            println(p*V[m,p,n,q])
            println(p*V[p,m,q,n])
            println("%%%%%%%%%%%")
        end
        return p*V[m,p,n,q]

    # If α excitation is two, it means m,n,p and q are all α.
    elseif αexcitation_level(D1, D2) == 2

        m,n = αexclusive_index(D1, D2)
        p,q = αexclusive_index(D2, D1)

        return p*(V[m,p,n,q] - V[m,q,n,p])

    # If α excitation is zero, it means m,n,p and q are all β.
    elseif αexcitation_level(D1, D2) == 0

        m,n = βexclusive_index(D1, D2)
        p,q = βexclusive_index(D2, D1)

        return p*(V[m,p,n,q] - V[m,q,n,p])
    end
end

end #module
