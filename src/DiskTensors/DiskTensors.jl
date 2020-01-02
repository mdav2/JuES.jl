module DiskTensors
include("DiskVectors.jl")
include("DiskMatrices.jl")
include("DiskFourTensors.jl")
export DiskVector
export blockfill!
export printdv
export dvwrite!
export dvread
export dvdot

export DiskMatrix
export dmread
export dmwrite
export dmdot

export DiskFourTensor
export d4read
#export dmelmult
end
