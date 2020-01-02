using HDF5

A = collect(reshape(1:120, 15, 8))
#h5write("/tmp/test.h5","mygroup2/A",A)
data = h5read("/tmp/test.h5","mygroup2/A",(2:3:15, 3:5)) #reading A[2:3:15,3:5] slice
println(data)

h5open("/tmp/test2.h5","w") do file
	A = rand(100,100)
	g = g_create(file,"mygroup")
	g["A","chunk",(5,5),"compress",3] = A
	attrs(g)["Description"] = "Just a single dataset"
end


