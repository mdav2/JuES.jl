# Input

The input file syntax in JuES is currently in a rudimentary place. The current style is,

```plain
molecule = """
<geometry>
"""
--
options = Dict(<psi_options>,
               <JuES_options)
--
<command>
```

Psi4 options must be wrapped with quotation marks `"scf_type" => "pk"` and JuES options are written with a preceeding
colon `:` indicating that they are Julia symbols. See `JuES.jl/test/` for some examples.

## Functions

```@docs
JuES.Input.run
```
```@docs
JuES.Input.read
```
```@docs
JuES.Input.exec
```
