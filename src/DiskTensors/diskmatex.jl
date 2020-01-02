using DiskTensors
using BenchmarkTools

N = 1000
a = DiskMatrix("/tmp/matrix1.h5",Float64,N,N,"w")
b = DiskMatrix("/tmp/matrix2.h5",Float64,N,N,"w")
blockfill!(a,1.0)
blockfill!(b,2.0)
println(@benchmark dmdot(a,b,201))
