module MatrixElement
"""
Module implementing matrix elements (currently just Hamiltonian) via
Slater's rules

Formulas and notation from Szabo and Ostlund p. 70

syntax --> <nparticle>Electron_<matrix>d<ndiff>
	   -->
	   --> e.g. for the one particle matrix element of a Hamiltonian betweeen
	   --> two determinants differing by 1 orbital:
	   --> 	OneElectron_Hd1
"""

using JuES.Wavefunction
using JuES.ConfigurationInteraction.DetOperations

export Hd0
export Hd1
export Hd2

function Hd0(a::Determinant, b::Determinant, ref::Wfn)
    """
    Σ <m|h|m> + 1/2 ΣΣ <mn||mn> 
    """
end
function Hd1(a::Determinant, b::Determinant, ref::Wfn)
    """
    differ m -> p
    <m|h|p> + Σ<mn||pn>
    """
end
function Hd2(a::Determinant, b::Determinant, ref::Wfn)
    """
    mn -> pq
    <mn||pq>
    """
end

end #module
