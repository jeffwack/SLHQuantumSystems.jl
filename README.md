# SLHQuantumSystems

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jeffwack.github.io/SLHQuantumSystems.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeffwack.github.io/SLHQuantumSystems.jl/dev/)
[![Build Status](https://github.com/jeffwack/SLHQuantumSystems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeffwack/SLHQuantumSystems.jl/actions/workflows/CI.yml?query=branch%3Amain)

SLHQuantumSystems is a Julia package for creating and combining open
quantum systems using the [SLH framework](https://arxiv.org/abs/1611.00375).

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

