using LinearAlgebra
using HDF5
struct DiskMatrix
	"""
	data structure for storing and accessing matrices on disk
	"""
	fname::String
	gname::String
	dname::String
	dtype::Union{Type{Float32},Type{Float64}}
	szx::Int #number of elements in x-axis
	szy::Int #number of elements in y-axis
end

function DiskMatrix(fname::String,dtype::Type,szx::Int,szy::Int,mode::String="r+")
	"Constructor for disk matrices"
	file = h5open(fname,mode)
	if (mode == "w") | (mode == "w+")
		dataset = d_create(file,"data",datatype(dtype),dataspace(szx))
	else
		dataset = file[fieldname]
	end
	DiskMatrix(fname,"data","data",dtype,szx,szy)
end

function blockfill!(dmat::DiskMatrix,val)
	"""
	Fill a DiskMatrix with a given value.
	TODO: replace with buffered I/O
	"""
	A = zeros(Float64,dmat.szx,dmat.szy)
	A .= val
	h5write(dmat.fname,"$dmat.dname",A)
end
function elmult()
end
function dmdot(dmat1::DiskMatrix,dmat2::DiskMatrix,buffsize)
	"""Buffered matrix matrix 'dot' product (accumulate multiply)"""
	chunks = cld(dmat1.szx,buffsize)
	tsum = 0.0
	for i in 1:1:dmat1.szx
		for chunk in 1:1:chunks-1
			tsum += dot(dmread(dmat1,(chunk-1)*buffsize+1:chunk*buffsize,i),
						dmread(dmat2,(chunk-1)*buffsize+1:chunk*buffsize,i))
		end
		tsum += dot(dmread(dmat1,(chunks-1)*buffsize+1:dmat1.szx,i),
					dmread(dmat2,(chunks-1)*buffsize+1:dmat1.szx,i))
	end
	return tsum
end
function dmdot(dmat1::DiskMatrix,dmat2::DiskMatrix)
	"""matrix matrix 'dot' product (accumulate multiply) - very slow!"""
	tsum = 0.0
	for i in 1:1:dmat1.szx
		for j in 1:1:dmat2.szy
			tsum += dmread(dmat1,i,j)[1]*dmread(dmat2,i,j)[1]
		end
	end
	return tsum
end
function dmread(dmat::DiskMatrix,posx,posy)
	fid = h5open(dmat.fname,"r")
	return fid["$dmat.dname"][posx,posy]
end
function dmwrite(dmat::DiskMatrix,posx,posy,val)
	fid = h5open(dmat.fname,"r+")
	fid["$dmat.dname"][posx,posy] = val
end
