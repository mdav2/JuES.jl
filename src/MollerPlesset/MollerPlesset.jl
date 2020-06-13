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
#using JuES.Direct
using JuES.IntegralTransformation
using JuES.DF
using JuES.Output
using JuES
using TensorOperations

export do_rmp2
export do_ump2
export do_direct_rmp2
export do_df_rmp2

function print_header()
    @output "================================================================================\n" 
    @output "|   Moller-Plesset Perturbation Theory                                         |\n"
    @output "|       module written by M.M. Davis                                           |\n"
    @output "================================================================================\n" 
end

include("RMP2.jl")
include("UMP2.jl")
#include("DirectRMP2.jl")
include("DF-RMP2.jl")

end #module
