# Coupled Cluster

```@docs
JuES.CoupledCluster
```

## Conventional Coupled Cluster 
### AutoRCCSD
Module that performs Restricted CCSD using auto factorized equations.
```@docs
JuES.CoupledCluster.AutoRCCSD
```
```@docs
JuES.CoupledCluster.AutoRCCSD.update_energy
```

### CCD
```@docs
JuES.CoupledCluster.RCCD
```
```@docs
JuES.CoupledCluster.RCCD.do_rccd
```

### CCSD
```@docs
JuES.CoupledCluster.RCCSD
```
```@docs
JuES.CoupledCluster.RCCSD.do_rccsd
```

## Density Fitted Coupled Cluster

### DFRCCD [BROKEN]
Module that performs density fitted CCD. This module is not operational, but should be fixed when the integrals
code is incorporated.

## GPU Accelerated Coupled Cluster

### GPRCCD
Module that performs GPU accelerated CCD. This module requires a functioning CUDA setup and an NVIDIA GPU.
