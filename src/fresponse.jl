#=
This file provides methods for evaluating the frequency domain response of a linear system.

For the input-to-output transfer matrix G(ω) as a plain array, use ControlSystems.jl's
`freqresp(sys, ωs)` directly — it works for QuantumStateSpace via the `ssdata` implementation
in abcd.jl.

`transfer_matrix` also includes the states, as outputs (the response of a mode to an input)
and as inputs (a unit excitation of a state, such as a force on a mirror), and is indexed by
name: Strings for ports and Symbols for states.
=#

using ControlSystems: freqresp, ss

"""
    TransferMatrix

Frequency response between the inputs, outputs and states of a `QuantumStateSpace`.

`G[:, :, k]` is the matrix at `freqs[k]`. Rows are the output ports followed by the states,
columns are the input ports followed by the states. Index by name, `T["out_p", "in_x"]`,
`T[:arm_x, "in_x"]`, or by integer; see [`transfer_matrix`](@ref).
"""
struct TransferMatrix
    freqs :: Vector{Float64}                    # angular frequencies [rad/s]
    G     :: Array{ComplexF64, 3}               # (nrows, ncols, nω)
    rows  :: Vector{Union{String, Symbol}}      # output ports, then states
    cols  :: Vector{Union{String, Symbol}}      # input ports, then states
end

"""
    transfer_matrix(sys::QuantumStateSpace, freqs) → TransferMatrix

Responses between every input, output and state of a numeric `sys` at the angular
frequencies `freqs` [rad/s]:

    [C; I] (iω - A)⁻¹ [B I] + [D 0; 0 0]

Ports are named by Strings and states by Symbols, each expanded to its two entries
(`"out_x"`, `"out_p"`, `:arm_x`, `:arm_p` in the quadrature basis). A name without its
suffix selects both entries:

    T = transfer_matrix(sys, ω)
    T["out_p", "in_x"]   # vector over frequency
    T[:arm_x, "in"]      # 1×2×nω array
"""
function transfer_matrix(sys::QuantumStateSpace, freqs::AbstractVector{<:Real})
    A, B, C, D = ssdata(sys)
    n, nin, nout = nstates(sys), ninputs(sys), noutputs(sys)
    augmented = ss(A, [B I(n)], [C; I(n)], [D zeros(nout, n); zeros(n, nin + n)])
    states = state_names(sys)
    return TransferMatrix(
        collect(Float64, freqs), freqresp(augmented, freqs),
        [port_names(sys.outputs, sys.basis); states],
        [port_names(sys.inputs, sys.basis); states]
    )
end

# The indices of `name` in `names`: an exact match, or else every entry that extends it by
# whole `_`-separated words, so "out" selects "out_x" and "out_p" but "ou" selects nothing
function lookup(names, name)
    i = findall(isequal(name), names)
    isempty(i) || return i
    i = findall(n -> n isa typeof(name) && startswith(string(n), string(name) * "_"), names)
    isempty(i) && error("no $(name isa String ? "port" : "state") named $(repr(name)); ports are Strings and states are Symbols. Available: $names")
    return i
end

# A vector over frequency for a single entry, otherwise an array
function select(M, rows, cols, i, j)
    ii, jj = lookup(rows, i), lookup(cols, j)
    return length(ii) == length(jj) == 1 ? M[only(ii), only(jj), :] : M[ii, jj, :]
end

const Name = Union{String, Symbol}

Base.getindex(T::TransferMatrix, i::Int, j::Int) = T.G[i, j, :]
Base.getindex(T::TransferMatrix, i::Name, j::Name) = select(T.G, T.rows, T.cols, i, j)

function Base.show(io::IO, T::TransferMatrix)
    f1 = round(T.freqs[1] / (2π), sigdigits = 3)
    f2 = round(T.freqs[end] / (2π), sigdigits = 3)
    println(io, "TransferMatrix: $(size(T.G, 1))×$(size(T.G, 2)) over $(length(T.freqs)) frequencies ($f1 – $f2 Hz)")
    println(io, "  To:   $(join(repr.(T.rows), ", "))")
    println(io, "  From: $(join(repr.(T.cols), ", "))")
end

function symbfresponse(sys::QuantumStateSpace)
    @variables s
    iden = Matrix{Int}(I, size(sys.A)...)
    G = inv(s*iden - sys.A)
    return simplify.(sys.C*G*sys.B + sys.D)
end
