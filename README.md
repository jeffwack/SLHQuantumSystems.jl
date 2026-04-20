# SLHQuantumSystems

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jeffwack.github.io/SLHQuantumSystems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeffwack.github.io/SLHQuantumSystems.jl/dev/)
[![Build Status](https://github.com/jeffwack/SLHQuantumSystems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeffwack/SLHQuantumSystems.jl/actions/workflows/CI.yml?query=branch%3Amain)

SLHQuantumSystems is a Julia package for creating and combining open
quantum systems using the [SLH framework](https://arxiv.org/abs/1611.00375).

## Scope

SLHQuantumSystems is a standalone package for:

- Creating SLH triples `(S, L, H)` with symbolic parameters
  ([Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)) and
  quantum operators
  ([SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl)).
- Composing SLH systems via `concatenate` and `feedbackreduce`.
- Converting **linear-bosonic** SLH systems to state-space form
  (`QuantumStateSpace <: ControlSystems.AbstractStateSpace`) for
  frequency-domain analysis with ControlSystems.jl (`freqresp`, Bode, etc.).
- Physical-units-aware output via mode metadata (`OpticalMode`,
  `MechanicalMode`) and SI conversion of spectral densities.

### Out of scope (by design)

- **Wiring-diagram primitives and port-based composition** (`Mirror`,
  `BeamSplitter`, `link!`, cavity detection, …) — these live in
  [Breadboard](https://github.com/jeffwack/Breadboard).
- **Non-bosonic Hilbert spaces** in the state-space conversion. Atomic
  (`NLevelSpace`) models compose at the SLH layer but are not linearized
  through the ABCD pipeline.
- **Hierarchical parameter storage.** SLH models use a flat parameter
  dict; names are promoted during `concatenate` to avoid collisions, but
  the parameter dict itself has no nested structure. Hierarchical / graph
  structure is the wiring layer's concern.
- **YAML/JSON model specs.** Model serialization belongs in Breadboard.
- **Time-domain simulation.** Use
  [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl),
  [QuantumCumulants.jl](https://github.com/qojulia/QuantumCumulants.jl),
  or [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl).

## Installation

This package is registered in the Julia registry, you can install it with 

```julia
|pkg> add SLHQuantumSystems
```

## Contributing
Contributions to this package are welcome! Before making a PR, you should:
- run the package tests
```julia
|pkg> test
```
- build the documentation
```julia
|pkg> activate ./docs/
|julia> using LiveServer; servedocs()
```

