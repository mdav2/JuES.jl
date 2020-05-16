# JuES : Documentation
## Overview and Motivations

JuES (pronounced "juice") is a programming environment for writing arbitrary electronic structure and quantum chemical computations in the Julia programming language. Julia shows a lot of promise as a language for scientific computing, with many fields writing domain-specific applications in Julia. This project is intended to demonstrate some ways of working in this language, and showcase a proposed style of programming for expansion into a complete set of electronic structure programs.

In the benchmark and test folders, there are some files that act as examples for using the program. 

## Installing JuES
These instructions are lifted from [this helpful site](https://tlienart.github.io/pub/julia/dev-pkg.html). 
In addition to the instructions below, there are some dependencies required. Please raise an issue if there is any difficulty with dependencies.
Make a directory where you will be placing JuES. I'll use `<DEVDIR>` to represent that directory. Clone JuES.

### Linking to Psi4Numpy
> **(!) Temporary** Due to a conflict between Numpy and MKL, if you want to use Intel MKL, which can boost performance for tensor contractions, currently you must obtain MKL.jl from [this repo](https://github.com/fgerick/MKL.jl) NOT the JuliaComputing master branch (currently). 

Please do the following to make the Psi4Numpy interface visible to Julia.
```
conda create -n p4env python=3.7 psi4 -c psi4
```
To get the path to the Python executable,
```
$ conda activate p4env
$ which python
```
Then make this python visible to Julia:
```julia-repl
julia>ENV["PYTHON"] = <path-to-p4env-python>
julia>] build PyCall
```

### Making Julia aware of JuES
```
mkdir <DEVDIR>
cd <DEVDIR>
git clone https://github.com/mdav2/JuES.jl.git
```
Now make the package manager Pkg aware of JuES' presence.
```
julia> ]
(v1.3) pkg> dev <DEVDIR>/JuES.jl
```
Now JuES should be visible to Julia! To test this, `cd` and run,
```
$ julia
julia> using JuES
```
You should see,
```
[ Info: Precompiling JuES [9237668d-08c8-4784-b8dd-383aa52fcf74]
```

I strongly encourage running the test suite immediately! Pkg once again makes this easy,
```
julia> ]
(v1.3) pkg> test JuES
```
or if you prefer
```
julia> import Pkg
julia> Pkg.test("JuES")
```
You will see a bunch of package versions spill out, and the REPL will hang while it runs tests. This takes about 2 minutes on a very slow machine (Intel Core m3 processor). Finally you will see an output from the test suite. If all is well, you will see something like,
```
Test Summary:             | Pass  Fail  Error  Total
JuES                      |   18                  18
```
The number of tests may be different, though. However, in the more likely case that something has an unfixed bug, in this case an issue with DiskTensors and an issue with MP2. "Fail" means the test ran without errors (syntax etc.) but did not produce the prescribed output. Usually this is a problem with the implementation of the method - please place an issue on GitHub! "Error" means the mistake was more boneheaded - usually because a breaking change was made without correctly modifying the test suite. 
```
Test Summary:             | Pass  Fail  Error  Total
JuES                      |   16     1      1     18
  Wavefunction            |    4                   4
  CISingles               |    1                   1
  DiskTensors             |    6            1      7
    Smoke                 |    3                   3
    Dot                   |    3            1      4
      DiskVector          |    1            1      2
      DiskMatrix          |    1                   1
      DiskFourTensor      |    1                   1
  CoupledCluster          |    1                   1
  Integral Transformation |    2                   2
  MP2                     |    2     1             3
```

