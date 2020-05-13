module PerturbativeTriples
using TensorOperations



function do_pT(T1 ::Array{Float64, 2}, T2::Array{Float64, 4}, Voovv::Array{Float64,4}, Vovvv::Array{Float64,4}, Vooov::Array{Float64,4}, fo::Array{Float64,1}, fv::Array{Float64,1}, chonky::Bool = false)

    if chonky

        # Build full resolvent
        Dd = [i + j + k - a - b - c for i = fo, j = fo, k = fo, a = fv, b = fv, c = fv]


        # Get connected terms
        t3c 
        

