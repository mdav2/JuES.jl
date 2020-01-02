using DiskTensors
a = DiskFourTensor("/tmp/dtens1.h5",Float64,10,10,10,10,"w")
blockfill!(a,1.0)
println(d4read(a,1,1,1,1))
