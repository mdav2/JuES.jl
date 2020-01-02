# JuES : Documentation
## Overview and Motivations
> **(!)** This project is preliminary. I'm still quite new to Julia, and any feedback (including recommending a complete refactoring!) is very welcome! 

JuES (pronounced like "juice") is a programming environment for writing arbitrary electronic structure and quantum chemical computations in the Julia programming language. Julia shows a lot of promise as a language for scientific computing, with many fields writing domain specific applications in Julia. This project is intended to demonstrate some ways of working in this language, and showcase a proposed style of programming for expansion into a complete set of electronic structure programs.
## Tests
Once installed, running tests is done like 
```
$ julia
julia > import Pkg
julia > Pkg.test("JuES")
```

## Design
### Integral backends
This project uses other electronic structure programs to compute basic quantities like the one- and two-electron Hamiltonian integrals, as well as their counterparts for other operators. There is currently no intention of writing an integrals code specifically for JuES.

As of version `alpha1` there is an interface to the Psi4 programs via the `psi4numpy` interface. This is simply the interface that I know the best, and should be extended in the future. There are plans for interfacing with the PySCF project and the NWChem project. If someone with knowledge wants to implement an interface to other programs e.g. Q-Chem, CFOUR, ORCA, Turbomole, or MOLPRO, they are very welcome. It is intended that the most robust interfaces should be to free and/or open source programs. 

### Naming
It is customary for Julia modules and data structures to have CamelCase names, such as `Wavefunction.jl`. Please follow this aesthetically pleasing convention! 

## Modules
This section contains a description of some preliminary modules in the JuES environment. Many aspects are aspirational, and the description is not so much a description of current functionality as a statement of intent.
### Wavefunction.jl
Wavefunction.jl is the foundational module of JuES. All JuES programs will make use of Wavefunction.jl at some point. This module is the point of interaction between integral backends (e.g. `psi4numpy`) and the JuES programming environment. An interface to an integral backend should produce a complete Wfn structure from the relevant sources. A Wfn structure is a representation of a reference determinant, such as a set of Hartree-Fock or DFT orbitals. 
> **(!) Help** This module is nearly feature complete, and has a skeleton test set. Fleshing out documentation and contributing tests would be a great help!
### Determinant.jl
General representation of a Slater determinant. Used in modules CISingles.jl. 
> **(!) Help** This module is severely lacking - excellent starting project for an experienced programmer.
### DiskTensors.jl
This module describes a way of storing vectors, matrices, and rank four tensors on disk in a convenient way for use in electronic structure computations. 
This module uses the HDF5 binary format for storing and accessing arrays. This was chosen for its convenient interface, support for compressed I/O, and good support in Julia. 
> **(!) Help** This module is about halfway feature-complete, but lacks a test suite and documentation.
### Davidson.jl
This module implements a simple Davidson solver, which currently has some unidentified bug. Use the IterativeSolvers.jl LOBPCG routine for in-core computations.
> **(!) Help** Rewriting the Davidson code is probably a good idea. A routine to collapse the trial vector subspace would make this module much more functional. A generalized implementation for non-symmetric matrices is required before EOM codes can be useful. 
### Direct.jl
This module contains necessary code for integral direct computations. Currently only interfaces to the Psi4 programs, but additions of interfaces are welcome. 
> **(!) Help** I don't yet understand how to code a reasonable integral direct program, so this will likely be neglected for some time without outside help.
### MatrixElement.jl
This module defines an interface for obtaining matrix elements for CI matrices.
>**(!) Help** This is just a skeleton at this point. Contributions to this module will greatly help a functioning FCI and arbitrary order CI code. 
### MollerPlesset.jl
Routines for Moller-Plesset perturbation theory computations are implemented here. Currently only an in-core RMP2 implementation exists. 
>**(!) Help** Disk based RMP2, as well as UMP2 codes would be an excellent contribution. Direct MP2 is also an excellent contribution.
### CISingles.jl
Specialized routines for computing configuration-interaction singles excited state wavefunctions are defined here. Corrections such as CIS(D) and variants defined here as well. Keep seperate from general CI code. Only in-core RCIS is implemented.
>**(!) Help** UCIS, disk-based, and direct implementations are excellent targets. 
### CoupledCluster.jl
Routines for computing ground state coupled cluster energies are contained here. Currently there is only RHF-CCD implemented. RHF-CCSD should be implemented soon, and UHF-CCD sometime after.
>**(!) Help** Refining the RHF-CCD implementation, or working on CCSD codes would be greatly appreciated. Once RHF-CCSD is complete, coding a perturbative triples correction would be a straightforward addition.
