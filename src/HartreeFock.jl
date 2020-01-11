module HartreeFock
using PyCall
const psi4 = PyNULL()
function __init__()
    copy!(psi4, pyimport("psi4"))
end

include("RHF.jl")
export RHF
export RHFWfn
end #module

