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

```@docs
JuES.Transformation.tei_transform
```

## DF.jl

```@docs
JuES.DF
```
