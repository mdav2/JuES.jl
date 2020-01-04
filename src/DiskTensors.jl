module DiskTensors
include("DiskVectors.jl")
include("DiskMatrices.jl")
include("DiskFourTensors.jl")
function ranger(inp::Union{UnitRange{Int64},Int64,Colon})
	if typeof(inp) == Int64
		return UnitRange(inp:inp)
	else
		return inp
	end
end
export ranger
export DiskVector
export blockfill!
export printdv
export dvwrite!
export dvread
export dvdot
#export Base.getindex
include("DiskVectorOps.jl")

export DiskMatrix
export dmread
export dmwrite
export dmdot

export DiskFourTensor
export d4read
#export dmelmult
end
