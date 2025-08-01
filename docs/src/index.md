# SLHQuantumSystems.jl

SLHQuantumSystems.jl is a Julia package for creating and composing open
quantum systems using the SLH framework. 

## Quick Start
Get started by running one of the examples!
```@repl
include("examples/cascadedcavities.jl")
```

## Overview of SLH systems

The SLH framework represents each open quantum systems with three components:
- **S**: Scattering matrix describing direct input-output coupling of external
  (bath) modes
- **L**: Coupling vector describing the interaction of the internal modes with
  the external modes 
- **H**: System Hamiltonian describing internal dynamics

### Component Library
The SLH framework enables you to create complicated quantum systems by combining
simple, reusable components
- Pre-built quantum components including:
  - Basic cavities
  - Squeezing cavities  
  - Radiation pressure cavities
  - Jaynes-Cummings QED cavity


## Dependencies

- [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) provides the symbolic algebra system for quantum operators

## References

- This package was inspired by [QNET](https://github.com/mabuchilab/QNET), a
  python package for working with SLH systems.

- [The SLH framework for modeling quantum input-output networks](https://arxiv.org/pdf/1611.00375)
