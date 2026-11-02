#=
This file provides the transfer function of a linear system, sampled over frequency or as a
symbolic function of s.

`transfer_function` includes the states, as outputs (the response of a mode to an input)
and as inputs (a unit excitation of a state, such as a force on a mirror), and is indexed by
name: Strings for ports and Symbols for states. For the input-to-output matrix G(ω) alone as
a plain array, ControlSystemsBase's `freqresp(sys, ωs)` also works on a numeric system.
=#

using ControlSystemsBase: freqresp, ss, ssdata, nstates, ninputs, noutputs

"""
    QuantumTransferFunction

Frequency response between the inputs, outputs and states of a `QuantumStateSpace`.

`G[:, :, k]` is the matrix at `freqs[k]`. Rows are the output ports followed by the states,
columns are the input ports followed by the states. Index by name, `T["out_p", "in_x"]`,
`T[:arm_x, "in_x"]`, or by integer; see [`transfer_function`](@ref).
"""
struct QuantumTransferFunction
    freqs :: Vector{Float64}                    # angular frequencies [rad/s]
    G     :: Array{ComplexF64, 3}               # (nrows, ncols, nω)
    rows  :: Vector{Union{String, Symbol}}      # output ports, then states
    cols  :: Vector{Union{String, Symbol}}      # input ports, then states
end

"""
    SymbolicQuantumTransferFunction

Transfer function between the inputs, outputs and states of a `QuantumStateSpace`, as a
function of the Laplace variable `s`. Rows and columns are named as in [`QuantumTransferFunction`](@ref).

Each entry is computed when it is indexed, `T["out_p", "in_x"]`, as a ratio of polynomials
in `s` that is not reduced to lowest terms; see [`transfer_function`](@ref).
"""
struct SymbolicQuantumTransferFunction
    s    :: Num
    M    :: Matrix{Num}                         # sI - A
    B    :: Matrix{Num}                         # [B I]
    C    :: Matrix{Num}                         # [C; I]
    D    :: Matrix{Num}                         # [D 0; 0 0]
    rows :: Vector{Union{String, Symbol}}       # output ports, then states
    cols :: Vector{Union{String, Symbol}}       # input ports, then states
    den  :: Base.RefValue{Num}                  # det(sI - A), computed on first use
end

"""
    transfer_function(sys::QuantumStateSpace, freqs) → QuantumTransferFunction

Responses between every input, output and state of a numeric `sys` at the angular
frequencies `freqs` [rad/s], that is at `s = iω`:

    [C; I] (iω - A)⁻¹ [B I] + [D 0; 0 0]

Ports are named by Strings and states by Symbols, each expanded to its two entries
(`"out_x"`, `"out_p"`, `:arm_x`, `:arm_p` in the quadrature basis). A name without its
suffix selects both entries:

    T = transfer_function(sys, ω)
    T["out_p", "in_x"]   # vector over frequency
    T[:arm_x, "in"]      # 1×2×nω array

Every parameter must have a value; `substitute` them first.
"""
function transfer_function(sys::QuantumStateSpace, freqs::AbstractVector{<:Real})
    A, B, C, D = numeric_ssdata(sys)
    n, nin, nout = nstates(sys), ninputs(sys), noutputs(sys)
    augmented = ss(A, [B I(n)], [C; I(n)], [D zeros(nout, n); zeros(n, nin + n)])
    states = state_names(sys)
    return QuantumTransferFunction(
        collect(Float64, freqs), freqresp(augmented, freqs),
        [output_names(sys); states],
        [input_names(sys); states]
    )
end

"""
    transfer_function(sys::QuantumStateSpace, s::Num) → SymbolicQuantumTransferFunction

The same responses as a symbolic function of the Laplace variable `s`, for a `sys` in the
quadrature basis that may be symbolic, numeric or a mix:

    @variables s
    T = transfer_function(sys, s)
    T["out_p", "in_x"]   # a ratio of polynomials in s

Each entry is computed by determinants when it is indexed, which is practical up to about
8 states. Substituting `s = iω` gives the frequency response at `ω`.
"""
function transfer_function(sys::QuantumStateSpace, s::Num)
    if !(sys.basis isa QuadratureBasis)
        return error("the symbolic transfer function is only implemented in the quadrature basis; convert with QuantumStateSpace(sys, QuadratureBasis())")
    end
    Symbolics.issym(Symbolics.unwrap(s)) || return error("s must be a single symbol, got $s")
    A, B, C, D = (realmatrix(M, sys.name) for M in ssdata(sys))
    if any(M -> any(x -> any(v -> isequal(v, s), Symbolics.get_variables(x)), M), (A, B, C, D))
        return error("$s is a parameter of $(sys.name); pass a different symbol")
    end
    n, nin, nout = nstates(sys), ninputs(sys), noutputs(sys)
    states = state_names(sys)
    return SymbolicQuantumTransferFunction(
        s, s * I - A,
        Matrix{Num}([B I(n)]), Matrix{Num}([C; I(n)]), Matrix{Num}([D zeros(nout, n); zeros(n, nin + n)]),
        [output_names(sys); states],
        [input_names(sys); states],
        Ref{Num}()
    )
end

# The ABCD matrices with every constant folded and a concrete eltype, or an error naming the
# parameters that remain
function numeric_ssdata(sys)
    mats = [substitute_matrix(M, Dict()) for M in ssdata(sys)]
    symbolic = [M for M in mats if !(eltype(M) <: Union{Float64, ComplexF64})]
    if !isempty(symbolic)
        params = unique(v for M in symbolic for x in M for v in Symbolics.get_variables(x))
        return error("$(sys.name) still has the parameters $(join(params, ", ")); substitute values for them first, or pass a symbol s for the symbolic transfer function")
    end
    return mats
end

# A matrix entry as a real Num. A physical system has real ABCD matrices in the quadrature basis,
# but entries can be stored as Complex{Num} or ComplexF64 with a zero imaginary part.
function realentry(x)
    x isa Complex && return issymzero(imag(x)) ? realentry(real(x)) : nothing
    x isa Real && !(x isa Num) && return Num(x)
    u = Symbolics.unwrap(x)
    return Symbolics.symtype(u) <: Real ? Num(u) : nothing
end

function realmatrix(M, name)
    return map(M) do x
        r = realentry(x)
        isnothing(r) && error("the symbolic transfer function needs real ABCD matrices, but $name has the complex entry $x")
        r
    end
end

# G = D + C(sI - A)⁻¹B entry by entry, using C_i (sI - A)⁻¹ B_j = -det([sI-A B_j; C_i 0]) / det(sI - A).
# Symbolics expands a determinant by minors, so there is no pivoting on symbolic entries.
function entry(T::SymbolicQuantumTransferFunction, i, j)
    if !isassigned(T.den)
        T.den[] = isempty(T.M) ? Num(1) : expand(det(T.M))
    end
    den = T.den[]
    border = [T.M T.B[:, j]; transpose(T.C[i, :]) 0]
    num = expand(T.D[i, j] * den - det(border))
    return issymzero(num) ? Num(0) : num / den
end

# The indices of `name` in `names`: an exact match, or else every entry that is `name` plus a
# basis suffix, so "out" selects "out_x" and "out_p" but not "out_aux_x", and "ou" selects nothing
function lookup(names, name)
    i = findall(isequal(name), names)
    isempty(i) || return i
    i = findall(n -> n isa typeof(name) && first(rsplit(string(n), "_"; limit = 2)) == string(name), names)
    isempty(i) && error("no $(name isa String ? "port" : "state") named $(repr(name)); ports are Strings and states are Symbols. Available: $names")
    return i
end

# A vector over frequency for a single entry, otherwise an array
function select(M, rows, cols, i, j)
    ii, jj = lookup(rows, i), lookup(cols, j)
    return length(ii) == length(jj) == 1 ? M[only(ii), only(jj), :] : M[ii, jj, :]
end

const Name = Union{String, Symbol}

Base.getindex(T::QuantumTransferFunction, i::Int, j::Int) = T.G[i, j, :]
Base.getindex(T::QuantumTransferFunction, i::Name, j::Name) = select(T.G, T.rows, T.cols, i, j)

# A single entry, otherwise a matrix of entries
Base.getindex(T::SymbolicQuantumTransferFunction, i::Int, j::Int) = entry(T, i, j)
function Base.getindex(T::SymbolicQuantumTransferFunction, i::Name, j::Name)
    ii, jj = lookup(T.rows, i), lookup(T.cols, j)
    return length(ii) == length(jj) == 1 ? entry(T, only(ii), only(jj)) : [entry(T, a, b) for a in ii, b in jj]
end

function Base.show(io::IO, T::QuantumTransferFunction)
    f1 = round(T.freqs[1] / (2π), sigdigits = 3)
    f2 = round(T.freqs[end] / (2π), sigdigits = 3)
    println(io, "QuantumTransferFunction: $(size(T.G, 1))×$(size(T.G, 2)) over $(length(T.freqs)) frequencies ($f1 – $f2 Hz)")
    println(io, "  To:   $(join(repr.(T.rows), ", "))")
    println(io, "  From: $(join(repr.(T.cols), ", "))")
end

function Base.show(io::IO, T::SymbolicQuantumTransferFunction)
    println(io, "SymbolicQuantumTransferFunction: $(length(T.rows))×$(length(T.cols)) in $(T.s)")
    println(io, "  To:   $(join(repr.(T.rows), ", "))")
    println(io, "  From: $(join(repr.(T.cols), ", "))")
end
