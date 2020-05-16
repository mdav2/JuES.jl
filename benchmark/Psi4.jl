__precompile__(false)
module Psi4
using PyCall
const psi4 = PyNULL()
function __init__()
    copy!(psi4, pyimport("psi4"))
end
#export psi4

end
