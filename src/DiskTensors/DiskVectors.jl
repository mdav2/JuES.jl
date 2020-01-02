using HDF5
using LinearAlgebra

function todisk(fname::String,A::Array{Float64})

end

struct DiskVector
	"data structure for holding file/group/dataset handles
	and other information required to use a disk based vector"
	fname::String #file name
	gname::String #group name
	dname::String #dataset name
	dtype::Union{Type{Float32},Type{Float64}}
	size::Int #number of elements to be allocated
	compress::Bool
	chonk::Bool 
	chonksz::Tuple
end

function DiskVector(fname::String,dtype::Type,
					size::Int,mode::String="r+",fieldname::String="data",
					compress::Bool=false,chonk::Bool=false,chonksz::Tuple=Tuple{Int,Int}((1,1)))
	"Constructor for DiskVector objects"
	file = h5open(fname,mode)
	#if opening w need to make group
	if (mode == "w") | (mode == "w+")
		dataset = d_create(file,fieldname,datatype(dtype),dataspace(size))
	else
		dataset = file[fieldname]
	end
	close(file)
	DiskVector(fname,fieldname,"data",dtype,size,compress,chonk,chonksz)
end

function blockfill!(dvec::DiskVector,val::Float64)
	"""
	Fill a DiskVector with a given value
	"""
	A = zeros(Float64,dvec.size)
	A .= val
	h5write(dvec.fname,"$dvec.dname",A)
end
function dvdot(dvec1::DiskVector,dvec2::DiskVector)
	"computes dot product of two DiskVectors - elementwise, very slow!"
	if (dvec1.size != dvec2.size)
		return false
	end
	tsum = 0.0
	for i in 1:1:dvec1.size
		tsum += dvread(dvec1,i)[1]*dvread(dvec2,i)[1]
	end
	return tsum
end
function dvdot(dvec1::DiskVector,dvec2::DiskVector,buffsize)
	"computes dot product of two DiskVectors - buffered, use this"
	chunks = cld(dvec1.size,buffsize)
	tsum = 0.0
	for chunk in 1:1:(chunks)
		tsum += dot(dvread(dvec1,(chunk-1)*buffsize+1:(chunk)*buffsize),
					dvread(dvec2,(chunk-1)*buffsize+1:(chunk)*buffsize))
				
	end
	return tsum
end
function dvwrite!(dvec::DiskVector,val::Float64,pos)
	"writes to a single position in a DiskVector"
	fid = h5open(dvec.fname,"r+")
	fid["$dvec.dname"][pos] = val
	close(fid)
end
function dvread(dvec::DiskVector,pos)
	"reads a single position in a DiskVector"
	fid = h5open(dvec.fname,"r")
	return fid["$dvec.dname"][pos]
end
function printdv(dvec)
	"print a DiskVector"
	println(h5read(dvec.fname,"$dvec.dname"))
end
