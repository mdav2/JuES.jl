"""
## IntegralTransformation
    JuES.IntegralTransformation

Module to handle integral transformations from AO to MO.

_Functions:_

get_eri  -> From a Wavefunction object, return a specified ERI array.

get_fock -> From a Wavefunction object, return the Fock matrix.

"""
module IntegralTransformation
using TensorOperations
using JuES.Wavefunction

export get_eri
export get_fock

"""
## get_eri 
    JuES.IntegralTransformation.get_eri(wfn, eri_string, notation="phys")

From a Wavefunction object, return a specified ERI array.

_Arguments:_

wfn        -> Wavefunction object.

eri_string -> String with length 4 identifying the type of ERI. 
              Characters must be (o, O, v, V). Each indicating Occupied and Virtual for ALPHA and beta.

notation   -> OPTIONAL. Values: "chem" or "phys". Return the array in Chemist's or Physicists' notation.

"""
function get_eri(wfn::Wfn, eri_string::String, notation::String = "phys")

    # Size of eri_string must be 4.
    if sizeof(eri_string) != 4
        error("Invalid string given to JuES.IntegralTransformation.get_eri: $eri_string")
    end

    # The computations below assumed Chemist's notation. Thus, the string is modified if Phys is requested
    if notation == "phys"
        eri_string = eri_string[[1,3,2,4]]
    end

    C = []
    # Get C1, C2, C3, C4 for the integral transformation
    for s in eri_string
        if s == 'o'
            push!(C, wfn.Cbo)
        elseif s == 'O'
            push!(C, wfn.Cao)
        elseif s == 'v'
            push!(C, wfn.Cbv)
        elseif s == 'V'
            push!(C, wfn.Cav)
        end
    end

    C1, C2, C3, C4 = C
    gao = wfn.ao_eri
    @tensoropt V[i,a,j,b] := C4[σ,b]*C3[λ,j]*C2[ν,a]*C1[μ,i]*gao[μ,ν,λ,σ]

    if notation == "phys"
        @tensor V[i,j,a,b] := V[i,a,j,b]
    end

    return V
end

"""
## get_fock
    JuES.IntegralTransformation.get_fock(wfn, spin ="alpha")

From a Wavefunction object, return the Fock matrix.

_Arguments:_

wfn  -> Wavefunction object.

spin -> OPTIONAL. String indicating the spin of the Fock matrix. Accept: "alpha", "a" or "up" for alpha arrays and "Beta", "b", "down" for beta arrays.
        Case insensitive. 

"""
function get_fock(wfn::Wfn, spin::String = "alpha")

    if lowercase(spin) in ["alpha", "up", "a"]
        C  = wfn.Ca
        Co = wfn.Cao
    elseif lowercase(spin) in ["beta", "down", "b"]
        C  = wfn.Cb
        Co = wfn.Cbo
    else
        error("Invalid Spin option given to JuES.IntegralTransformation.get_fock: $spin")
    end

    gao = wfn.ao_eri
    hao = wfn.hao

    @tensoropt (p=>100x, q=>100x, k=>x, μ=>100x, ν=>100x, λ=>100x, σ=>100x) begin
        f[p,q] := C[μ,p]*C[ν,q]*hao[μ,ν] 
        f[p,q] += 2*C[μ,p]*C[ν,q]*Co[λ,k]*Co[σ,k]*gao[μ,ν,λ,σ]
        f[p,q] -= C[μ,p]*C[ν,q]*Co[λ,k]*Co[σ,k]*gao[μ,λ,ν,σ]
    end

    return f
end

end # Module
