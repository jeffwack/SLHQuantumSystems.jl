#=
This file provides methods for evaluating the frequency domain response of a linear system.

For the full input-to-output transfer matrix G(ω), use ControlSystems.jl's `freqresp(sys, ωs)`
directly — it works for QuantumStateSpace via the `ssdata` implementation in abcd.jl.

`fresponse_state2output` handles the quantum-specific case of state-to-output transfer:
the response C(iω - A)⁻¹eⱼ from a unit excitation on state quadrature j to output quadrature k.
This has no ControlSystems.jl equivalent (state excitations are not inputs). It is implemented
by constructing a dummy system with B = eⱼ and D = 0 and delegating to `freqresp`.
=#

using ControlSystems: freqresp

"""
    fresponse_state2output(sys::QuantumStateSpace, freqs, from::Int, to::Int) → Vector{ComplexF64}

Transfer from state quadrature index `from` to output quadrature index `to` over `freqs` [rad/s].

Computes C[to,:] * (iω·I - A)⁻¹ * eₓ for each ω, where eₓ is a unit vector selecting
state `from`. For QuadratureBasis systems, states are ordered (x₁, p₁, x₂, p₂, …).

Implemented by passing a dummy system with B = eₓ, D = 0 to `freqresp`, so it inherits
all of ControlSystems.jl's numerics (Hessenberg form, etc.) rather than using a naive
matrix inverse.
"""
function fresponse_state2output(sys::QuantumStateSpace, freqs::AbstractVector{<:Real},
                                from::Int, to::Int)
    n = nstates(sys)
    B_pick = zeros(ComplexF64, n, 1)
    B_pick[from, 1] = 1
    D_zero = zeros(ComplexF64, noutputs(sys), 1)
    dummy = QuantumStateSpace(sys.name, sys.subspaces, sys.parameters,
                              ["state_$from"], sys.outputs,
                              sys.A, B_pick, sys.C, D_zero,
                              sys.timeevol, sys.basis)
    G = freqresp(dummy, freqs)   # (nout, 1, nω)
    return G[to, 1, :]
end

function symbfresponse(sys::QuantumStateSpace)
    @variables s
    iden = Matrix{Int}(I, size(sys.A)...)
    G = inv(s*iden - sys.A)
    return simplify.(sys.C*G*sys.B + sys.D)
end

