# Backend

## Wavefunction.jl
`Wavefunction.jl` defines the `Wfn` type and constructors. These
constructors take in a Psi4 wavefunction object (PyObject) and convert
necessary information to Julia data structures.
```@docs
JuES.Wavefunction
```

```@docs
JuES.Wavefunction.Wfn(wfn)
```

```@docs
JuES.Wavefunction.Wfn
```
## Transformation.jl

```@docs
JuES.Transformation
```
## IntegralTransformation.jl
This module contain functions used to convert AO integrals into MO integrals.
It can be used to get ERI arrays and Fock matrices.
```@docs
JuES.IntegralTransformation
```

```@docs
JuES.IntegralTransformation.get_eri
```

```@docs
JuES.IntegralTransformation.get_fock
```

```@docs
JuES.Transformation.tei_transform
```

## DF.jl

```@docs
JuES.DF
```
