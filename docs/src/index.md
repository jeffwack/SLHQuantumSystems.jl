# SLHQuantumSystems.jl

SLHQuantumSystems.jl is a Julia package for creating and composing open
quantum systems using the SLH framework. 

## Quick Start
Get started by running one of the examples!
```@repl
include("examples/cascadedcavities.jl")
```

## Overview of SLH systems

The SLH framework represents each open quantum systems by a triple containing:
- **S**: Scattering matrix describing direct input-output coupling of external
  (bath) modes
- **L**: Coupling vector describing the interaction of internal modes with
  external modes 
- **H**: Hamiltonian describing internal dynamics


## Scope of this package

SLHQuantumSystems.jl is for:

- Creating SLH triples `(S, L, H)` with symbolic parameters (Symbolics.jl)
  and quantum operators (SecondQuantizedAlgebra.jl).
- Composing named SLH 'blocks' via `concatenate` and `feedbackreduce`. Input,
  output, operator, and parameter names are promoted during composition.
- Converting linear-bosonic SLH systems to state-space form
  `QuantumStateSpace <: ControlSystems.AbstractStateSpace`.

### Component Library

A small set of reusable components ([`cavity`](@ref), [`squeezing_cavity`](@ref)). 

## Dependencies

- [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) provides the symbolic algebra system for quantum operators

## References

- This package was inspired by [QNET](https://github.com/mabuchilab/QNET), a
  python package for working with SLH systems.

- [The SLH framework for modeling quantum input-output networks](https://arxiv.org/pdf/1611.00375) [Combes_Kerckhoff_Sarovar_2017](@cite)

```@bibliography
```
