# System Image Generation

If you are not developing or frequently updating JuES, we recommend compiling the package into a custom system image
to (greatly) reduce startup times. This is straightforward with the PackageCompiler.jl package.

Navigate to `JuES.jl/src/SysImg` and execute the following instruction.
```julia
using PackageCompiler
create_sysimage(:JuES,sysimage_path="sys_JuES.so",precompile_execution_file="executor.jl")
```

After working for a couple of minutes, JuES and it's dependencies will have been collected and compiled into the system image
located in `JuES.jl/src/SysImg/sys_JuES.so`. With this, the start times for JuES can be massively reduced. Simply call Julia with
the following command to use this system image.

```shell
$ julia --sysimage <path-to-JuES>/src/SysImg
```

We observe about a 10X reduction in load times with this method. 

