using HDF5
using LinearAlgebra
import Base.getindex
import Base.setindex!

"data structure for holding file/group/dataset handles
and other information required to use a disk based vector"
struct DiskVector
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
function DiskVector(
    fname::String,
    dtype::Type,
    size::Int,
    mode::String = "r+",
    fieldname::String = "data",
    compress::Bool = false,
    chonk::Bool = false,
    chonksz::Tuple = Tuple{Int,Int}((1, 1)),
)
    "Constructor for DiskVector objects"
    h5open(fname, mode) do file
        #if opening w need to make group
        if (mode == "w") | (mode == "w+")
            dataset = d_create(file, fieldname, datatype(dtype), dataspace(size))
        else
            dataset = file[fieldname]
        end
    end
    DiskVector(fname, fieldname, "data", dtype, size, compress, chonk, chonksz)
end
function getindex(dvec::DiskVector, i1::Int64)
    "defines A[1] accessor"
    h5open(dvec.fname, "r") do fid
        fid["$dvec.dname"][i1][1]
        #return out
    end
    #dvread(dvec,i1)[1]
end
function getindex(dvec::DiskVector, i1::UnitRange{Int64})
    "defines A[1:2] accessor"
    h5open(dvec.fname, "r") do fid
        fid["$dvec.dname"][i1]
        #return out
    end
    #dvread(dvec,i1)
end
function setindex!(dvec::DiskVector, val, ::Colon)
    slice = 1:dvec.size
    h5open(dvec.fname, "r+") do fid
        fid["$dvec.dname"][slice] = val
    end
end
function setindex!(dvec::DiskVector, val, i1::Int64)
    "defines A[1] = 1.0 syntax"
    h5open(dvec.fname, "r+") do fid
        fid["$dvec.dname"][i1] = val
    end
end
function setindex!(dvec::DiskVector, val, i1::UnitRange{Int64})
    "defines A[1:2] = [1.0,2.0] syntax"
    h5open(dvec.fname, "r+") do fid
        fid["$dvec.dname"][i1] = val
    end
end

"""
Fill a DiskVector with a given value
"""
function blockfill!(dvec::DiskVector, val::Float64)
    A = zeros(Float64, dvec.size)
    A .= val
    h5write(dvec.fname, "$dvec.dname", A)
end
