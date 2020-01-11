"""
    JuES.MollerPlesset

module for running MP2 energies on restricted and unrestricted HF references.
## methods
    do_rmp2 -> see docstring ?JuES.MollerPlesset.do_rmp2
    do_ump2 -> see docstring ?JuES.MollerPlesset.do_ump2
"""
module MollerPlesset

using JuES.Wavefunction
using JuES.DiskTensors

export do_rmp2
export do_ump2

include("RMP2.jl")
include("UMP2.jl")

end #module
