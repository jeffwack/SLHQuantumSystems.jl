"""
    SpectralDensityMatrix

Frequency-resolved noise spectral density matrix S(ω) of the outputs and states.

Stores S as a 3D array with shape (n, n, nω), following the same
axis convention as `freqresp`: `S[:, :, k]` is the full matrix at
`freqs[k]`.

Field `names` contains the output port names followed by the state names, indexed
like a [`TransferMatrix`](@ref): `sd["out_p", "out_p"]`, `sd[:arm_x, :arm_x]`.
"""
struct SpectralDensityMatrix
    freqs :: Vector{Float64}                    # angular frequencies [rad/s]
    S     :: Array{ComplexF64, 3}               # (n, n, nω)
    names :: Vector{Union{String, Symbol}}      # output ports, then states
end

# Integer indexing: returns vector over frequency
Base.getindex(sd::SpectralDensityMatrix, i::Int, j::Int) = sd.S[i, j, :]
Base.getindex(sd::SpectralDensityMatrix, i::Name, j::Name) = select(sd.S, sd.names, sd.names, i, j)

function Base.show(io::IO, sd::SpectralDensityMatrix)
    n  = size(sd.S, 1)
    nω = length(sd.freqs)
    f1 = round(sd.freqs[1]   / (2π), sigdigits=3)
    f2 = round(sd.freqs[end] / (2π), sigdigits=3)
    println(io, "SpectralDensityMatrix: $(n)×$(n) over $(nω) frequencies ($f1 – $f2 Hz)")
    println(io, "  Channels: $(join(repr.(sd.names), ", "))")
end

"""
    vacuum_noise(sys::QuantumStateSpace) → Matrix

Return the vacuum input noise covariance matrix (1/2)·I, sized for
the number of input channels (2 quadratures per port).
"""
function vacuum_noise(sys::QuantumStateSpace)
    nin = 2 * length(sys.inputs)
    return Matrix{Float64}(I, nin, nin) / 2
end

"""
    spectral_density(sys::QuantumStateSpace, freqs; S_in=vacuum_noise(sys))
    → SpectralDensityMatrix

Compute the noise spectral density matrix of the outputs and states

    S(ω) = G(ω) · S_in · G†(ω)

where G(ω) is the response of the outputs and states to the input ports, from
[`transfer_matrix`](@ref).

# Arguments
- `sys`:   A numeric (parameter-substituted) QuantumStateSpace.
- `freqs`: Vector of angular frequencies [rad/s].
- `S_in`:  Input noise covariance. Defaults to vacuum (1/2)·I.
"""
function spectral_density(sys::QuantumStateSpace, freqs::AbstractVector{<:Real};
                           S_in = vacuum_noise(sys))
    T = transfer_matrix(sys, freqs)
    G = T.G[:, 1:ninputs(sys), :]      # (nout + nstates, nin, nω)
    n, _, nω = size(G)

    S = Array{ComplexF64}(undef, n, n, nω)
    for k in 1:nω
        Gk = G[:, :, k]
        S[:, :, k] = Gk * S_in * adjoint(Gk)
    end

    return SpectralDensityMatrix(T.freqs, S, T.rows)
end

