using HDF5

println(h5read("/tmp/testdv.h5","data/data"))
#fid=h5open("/tmp/testdv.h5","r")
#g=fid["data"]
#dset=g["data"]
#println(dset[2,2])
#close(fid)
