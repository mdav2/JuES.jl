using Davidson
using LinearAlgebra
using BenchmarkTools

n = 100
A = zeros(n, n)
A += I(n)
for i in UnitRange(1, n)
    A[i, i] = i
end
A += randn(n, n) * 0.05

A += transpose(A)
A ./ 2.0
A = Symmetric(A)

println(@benchmark eigdav(A, 3, 4, 40, 1 * 10^-6))
#println(eigdav(A,3,4,40,1*10^-6))
#val,vec = eigen(A)
#print(val[1:2])
