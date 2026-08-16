# SLHQuantumSystems

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jeffwack.github.io/SLHQuantumSystems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeffwack.github.io/SLHQuantumSystems.jl/dev/)
[![Build Status](https://github.com/jeffwack/SLHQuantumSystems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeffwack/SLHQuantumSystems.jl/actions/workflows/CI.yml?query=branch%3Amain)

SLHQuantumSystems is a Julia package for creating and composing open
quantum systems using the [SLH framework](https://arxiv.org/abs/1611.00375).

## Scope

SLHQuantumSystems is for:

- Creating SLH triples `(S, L, H)` with symbolic parameters
  ([Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)) and
  quantum operators
  ([SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl)).
- Composing SLH systems via `concatenate` and `feedbackreduce`.
- Converting **linear-bosonic** SLH systems to state-space form
  (`QuantumStateSpace <: ControlSystems.AbstractStateSpace`) for
  frequency-domain analysis with ControlSystems.jl (`freqresp`, Bode, etc.).

## Installation

```julia
|pkg> add SLHQuantumSystems
```

## Contributing
Contributions are welcome! Before making a PR, you should:
- run the package tests
```julia
|pkg> test
```
- build the documentation
```julia
|pkg> activate ./docs/
|julia> using LiveServer; servedocs(skip_dir="docs/src/generated")
```

