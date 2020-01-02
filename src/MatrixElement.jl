module MatrixElement
"""Module implementing matrix elements (currently just Hamiltonian) via
Slater's rules
syntax --> <nparticle>Electron_<matrix>d<ndiff>
	   -->
	   --> e.g. for the one particle matrix element of a Hamiltonian betweeen
	   --> two determinants differing by 1 orbital:
	   --> 	OneElectron_Hd1
"""

using JuES.Determinant
using JuES.Wavefunction

function OneElectron_Hd0(a::SlaterDeterminant,b::SlaterDeterminant,ref::Wfn)
end
function TwoElectron_Hd0(a::SlaterDeterminant,b::SlaterDeterminant,ref::Wfn)
end
function OneElectron_Hd1(a::SlaterDeterminant,b::SlaterDeterminant,ref::Wfn)
end
function TwoElectron_Hd1(a::SlaterDeterminant,b::SlaterDeterminant,ref::Wfn)
end
function OneElectron_Hd2(a::SlaterDeterminant,b::SlaterDeterminant,ref::Wfn)
end
function TwoElectron_Hd2(a::SlaterDeterminant,b::SlaterDeterminant,ref::Wfn)
end

end #module
