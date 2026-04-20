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


## Scope of this package

SLHQuantumSystems.jl is a standalone package for:

- Creating SLH triples `(S, L, H)` with symbolic parameters (Symbolics.jl)
  and quantum operators (SecondQuantizedAlgebra.jl).
- Composing named SLH 'blocks' via `concatenate` and `feedbackreduce`. Input,
  output, operator, and parameter names are promoted hierarchically during
  composition to avoid collisions.
- Converting **linear-bosonic** SLH systems to state-space form
  (`QuantumStateSpace <: ControlSystems.AbstractStateSpace`) for
  frequency-domain analysis with ControlSystems.jl (`freqresp`, Bode, …).
- Physical-units-aware output via mode metadata (`OpticalMode`,
  `MechanicalMode`) and SI conversion of spectral densities.

### Out of scope (by design)

- **Wiring-diagram primitives and port-based composition** (`Mirror`,
  `BeamSplitter`, `link!`, cavity detection). These live in Breadboard.
- **Non-bosonic Hilbert spaces** in the state-space conversion. Atomic
  (`NLevelSpace`) models compose at the SLH layer but are not linearized
  through the ABCD pipeline.
- **Hierarchical parameter storage.** SLH models carry a flat parameter
  dict. Promotion during `concatenate` keeps names unique, but the dict
  has no nested structure; hierarchy is the wiring layer's concern.
- **YAML/JSON model specs** and **time-domain simulation**.

### Component Library

A small set of reusable components ([`cavity`](@ref),
[`squeezing_cavity`](@ref)). The library is intentionally thin; richer
physical primitives are planned in the Breadboard wiring-diagram layer.


## Dependencies

- [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) provides the symbolic algebra system for quantum operators

## References

- This package was inspired by [QNET](https://github.com/mabuchilab/QNET), a
  python package for working with SLH systems.

- [The SLH framework for modeling quantum input-output networks](https://arxiv.org/pdf/1611.00375) [Combes_Kerckhoff_Sarovar_2017](@cite)

```@bibliography
```
