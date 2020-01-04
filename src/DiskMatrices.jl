using LinearAlgebra
using HDF5
import Base.getindex
import Base.setindex!
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
# >>> overload getindex -> A[i,j] for DiskMatrices
function getindex(dmat::DiskMatrix,
				  i1::Union{UnitRange{Int64},Colon},
				  i2::Union{UnitRange{Int64},Colon})
	h5open(dmat.fname, "r") do fid
		fid["$dmat.dname"][i1,i2]
	end
end
function getindex(dmat::DiskMatrix,
				  i1::Union{UnitRange{Int64},Int64,Colon},
				  i2::Union{UnitRange{Int64},Int64,Colon})
	getindex(dmat,ranger(i1),ranger(i2))
end
# <<< 
#
# >>> overload setindex! -> A[i,j] = val for DiskMatrices
function setindex!(dmat::DiskMatrix,val,i1::Union{UnitRange{Int64},Colon},i2::Union{UnitRange{Int64},Colon})
	h5open(dmat.fname, "r+") do fid
		fid["$dmat.dname"][i1,i2] = val
	end
end
function setindex!(dmat::DiskMatrix,val,i1::Int64,i2::UnitRange{Int64})
	setindex!(dmat,val,UnitRange(i1:i1),i2)
end
function setindex!(dmat::DiskMatrix,val,i1::UnitRange{Int64},i2::Int64)
	setindex!(dmat,val,i1,UnitRange(i2:i2))
end
function setindex!(dmat::DiskMatrix,val,i1::Int64,i2::Int64)
	setindex!(dmat,val,UnitRange(i1:i1),UnitRange(i2:i2))
end
# <<<
function blockfill!(dmat::DiskMatrix,val)
	"""
	Fill a DiskMatrix with a given value.
	TODO: replace with buffered I/O
	"""
	A = zeros(Float64,dmat.szx,dmat.szy)
	A .= val
	h5write(dmat.fname,"$dmat.dname",A)
end
function dmdot(dmat1::DiskMatrix,dmat2::DiskMatrix,buffsize)
	"""Buffered matrix matrix 'dot' product (accumulate multiply)"""
	chunks = cld(dmat1.szx,buffsize)
	tsum = 0.0
	for i in 1:1:dmat1.szx
		for chunk in 1:1:chunks-1
			tsum += dot(dmat1[(chunk-1)*buffsize+1:chunk*buffsize,i],
						dmat2[(chunk-1)*buffsize+1:chunk*buffsize,i])
		end
		tsum += dot(dmat1[(chunks-1)*buffsize+1:dmat1.szx,i],
				    dmat2[(chunks-1)*buffsize+1:dmat1.szx,i])
	end
	return tsum
end
function dmdot(dmat1::DiskMatrix,dmat2::DiskMatrix)
	"""matrix matrix 'dot' product (accumulate multiply) - very slow!"""
	tsum = 0.0
	for i in 1:1:dmat1.szx
		for j in 1:1:dmat2.szy
			tsum += dmat1[i,j]*dmat2[i,j]
		end
	end
	return tsum
end
#function dmread(dmat::DiskMatrix,posx,posy)
#	h5open(dmat.fname, "r") do fid
#		fid["$dmat.dname"][posx,posy]
#	end
#end
#function dmwrite(dmat::DiskMatrix,posx,posy,val)
#	h5open(dmat.fname, "r+") do fid
#		fid["$dmat.dname"][posx,posy] = val
#	end
#end
