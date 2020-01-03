using HDF5
using LinearAlgebra
import Base.getindex
import Base.setindex!

function todisk(fname::String,A::Array{Float64})

end

struct DiskVector
	"data structure for holding file/group/dataset handles
	and other information required to use a disk based vector"
	#file
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
	h5open(fname,mode) do file
		#if opening w need to make group
		if (mode == "w") | (mode == "w+")
			dataset = d_create(file,fieldname,datatype(dtype),dataspace(size))
		else
			dataset = file[fieldname]
		end
	end
	DiskVector(fname,fieldname,"data",dtype,size,compress,chonk,chonksz)
end
function getindex(dvec::DiskVector, i1::Int64)
	"defines A[1] accessor"
	h5open(dvec.fname,"r") do fid
		fid["$dvec.dname"][i1][1]
		#return out
	end
	#dvread(dvec,i1)[1]
end
function getindex(dvec::DiskVector, i1::UnitRange{Int64})
	"defines A[1:2] accessor"
	h5open(dvec.fname,"r") do fid
		fid["$dvec.dname"][i1]
		#return out
	end
	#dvread(dvec,i1)
end
function setindex!(dvec::DiskVector, val,::Colon)
	slice = 1:dvec.size
	h5open(dvec.fname,"r+") do fid
		fid["$dvec.dname"][slice] = val
	end
end
function setindex!(dvec::DiskVector, val, i1::Int64)
	"defines A[1] = 1.0 syntax"
	h5open(dvec.fname,"r+") do fid
		fid["$dvec.dname"][i1] = val
	end
end
function setindex!(dvec::DiskVector, val, i1::UnitRange{Int64})
	"defines A[1:2] = [1.0,2.0] syntax"
	h5open(dvec.fname,"r+") do fid
		fid["$dvec.dname"][i1] = val
	end
end

function blockfill!(dvec::DiskVector,val::Float64)
	"""
	Fill a DiskVector with a given value
	"""
	A = zeros(Float64,dvec.size)
	A .= val
	h5write(dvec.fname,"$dvec.dname",A)
end
#function dvdot(dvec1::DiskVector,dvec2::DiskVector)
#	"computes dot product of two DiskVectors - elementwise, very slow!"
#	if (dvec1.size != dvec2.size)
#		return false
#	end
#	tsum = 0.0
#	for i in 1:1:dvec1.size
#		tsum += getindex(dvec1,i)[1]*getindex(dvec2,i)[1]
#		#tsum += dvec1[i][1]*dvec2[i][1]
#		#tsum += dvread(dvec1,i)[1]*dvread(dvec2,i)[1]
#	end
#	return tsum
#end
#function dvdot(dvec1::DiskVector,dvec2::DiskVector,buffsize)
#	"computes dot product of two DiskVectors - buffered, use this"
#	chunks = cld(dvec1.size,buffsize)
#	tsum = 0.0
#	for chunk in 1:1:(chunks)
#		tsum += dot(dvread(dvec1,(chunk-1)*buffsize+1:(chunk)*buffsize),
#					dvread(dvec2,(chunk-1)*buffsize+1:(chunk)*buffsize))
#				
#	end
#	return tsum
#end
function dvwrite!(dvec::DiskVector,val::Float64,pos)
	"writes to a single position in a DiskVector"
	h5open(dvec.fname,"r+") do fid
		fid["$dvec.dname"][pos] = val
	end
end
function dvread(dvec::DiskVector,pos)
	"reads a single position in a DiskVector"
	h5open(dvec.fname,"r") do fid
		fid["$dvec.dname"][pos]
		#return out
	end
end
function printdv(dvec)
	"print a DiskVector"
	println(h5read(dvec.fname,"$dvec.dname"))
end
