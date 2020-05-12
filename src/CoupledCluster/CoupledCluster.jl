"""
Module for running CC computations in Julia.

Implemented --> RCCD, RCCSD, DF-RCCD

"""
module CoupledCluster

using JuES.Wavefunction
using JuES.Transformation
using JuES.Output
using Printf
using Base.Threads
#using SharedArrays
#using Distributed
using TensorOperations
using LinearAlgebra
using Dates

function print_header()
    @output repeat("=",80)*"\n"
    @output "|   {:<74} |\n" "Coupled Cluster"
    @output "|       {:<70} |\n" "Module by M.M. Davis and G.J.R. Aroeira"
    @output repeat("=",80)*"\n"
end

include("RCCD.jl")
include("DF-RCCD.jl")
include("ROCCD.jl")
include("RCCSD.jl")
include("UCCSD.jl")
include("unf-RCCSD.jl")
end #module CC
