using Determinant

i_a = BitArray(undef,4)
i_b = BitArray(undef,4)
j_a = BitArray(undef,4)
j_b = BitArray(undef,4)

i_a[:] = [1,1,0,0]
i_b[:] = [1,0,0,1]
j_a[:] = [1,0,1,0]
j_b[:] = [1,1,0,0]

i = SlaterDeterminant(i_a,i_b)
j = SlaterDeterminant(j_a,j_b)
println(norbdiff(i,j))
println(orbdiff(i,j))
