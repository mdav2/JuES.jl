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

using JuES.Determinant
using JuES.Wavefunction

function OneElectron_Hd0(a::SlaterDeterminant, b::SlaterDeterminant, ref::Wfn)
    """
    Σ <m|h|m>
    """
end
function TwoElectron_Hd0(a::SlaterDeterminant, b::SlaterDeterminant, ref::Wfn)
    """
    1/2 ΣΣ <mn||mn>
    """
end
function OneElectron_Hd1(a::SlaterDeterminant, b::SlaterDeterminant, ref::Wfn)
    """
    <m|h|p>
    """
end
function TwoElectron_Hd1(a::SlaterDeterminant, b::SlaterDeterminant, ref::Wfn)
    """
    differ m -> p
    Σ<mn||pn>
    """
end
function OneElectron_Hd2(a::SlaterDeterminant, b::SlaterDeterminant, ref::Wfn)
    """
    0
    """
end
function TwoElectron_Hd2(a::SlaterDeterminant, b::SlaterDeterminant, ref::Wfn)
    """
    mn -> pq
    <mn||pq>
    """
end

end #module