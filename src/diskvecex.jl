using DiskTensors
using BenchmarkTools
using LinearAlgebra

N = 2000000
x = DiskVector("/tmp/testdv1.h5",Float64,N,"w")
y = DiskVector("/tmp/testdv2.h5",Float64,N,"w")
blockfill!(x,1.0)
blockfill!(y,1.0)
#println(@benchmark dvdot(x,y,100))
println(@benchmark dvdot(x,y,5000))
println(@benchmark dot(ones(N),ones(N)))
#println(dvread(x,1:4))
