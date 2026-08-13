using ControlSystems: freqresp

"""
    SpectralDensityMatrix

Frequency-resolved output noise spectral density matrix S_out(ω).

Stores S as a 3D array with shape (nout, nout, nω), following the same
axis convention as `freqresp`: `S[:, :, k]` is the full matrix at
`freqs[k]`.

Field `names` contains the quadrature-expanded output port names.
"""
struct SpectralDensityMatrix
    freqs :: Vector{Float64}          # angular frequencies [rad/s]
    S     :: Array{ComplexF64, 3}     # (nout, nout, nω)
    names :: Vector{String}           # quadrature-expanded output port names
end

# Integer indexing: returns vector over frequency
Base.getindex(sd::SpectralDensityMatrix, i::Int, j::Int) = sd.S[i, j, :]

# Named indexing: looks up port names
function Base.getindex(sd::SpectralDensityMatrix, i::String, j::String)
    ii = findfirst(==(i), sd.names)
    jj = findfirst(==(j), sd.names)
    isnothing(ii) && error("Port \"$i\" not found. Available: $(sd.names)")
    isnothing(jj) && error("Port \"$j\" not found. Available: $(sd.names)")
    return sd.S[ii, jj, :]
end

Base.getindex(sd::SpectralDensityMatrix, i::Symbol, j::Symbol) =
    sd[string(i), string(j)]

function Base.show(io::IO, sd::SpectralDensityMatrix)
    n  = size(sd.S, 1)
    nω = length(sd.freqs)
    f1 = round(sd.freqs[1]   / (2π), sigdigits=3)
    f2 = round(sd.freqs[end] / (2π), sigdigits=3)
    println(io, "SpectralDensityMatrix: $(n)×$(n) over $(nω) frequencies ($f1 – $f2 Hz)")
    println(io, "  Channels: $(join(sd.names, ", "))")
end

# Private: expand port names to quadrature pairs
function _quadrature_output_names(sys::QuantumStateSpace{TE, QuadratureBasis}) where TE
    names = String[]
    for port in sys.outputs
        push!(names, port * "_x")
        push!(names, port * "_p")
    end
    return names
end

function _quadrature_output_names(sys::QuantumStateSpace{TE, LadderBasis}) where TE
    names = String[]
    for port in sys.outputs
        push!(names, port * "_a")
        push!(names, port * "_adag")
    end
    return names
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

Compute the output noise spectral density matrix

    S_out(ω) = G(ω) · S_in · G†(ω)

where G(ω) = freqresp(sys, freqs) is the transfer matrix.

# Arguments
- `sys`:   A numeric (parameter-substituted) QuantumStateSpace.
- `freqs`: Vector of angular frequencies [rad/s].
- `S_in`:  Input noise covariance. Defaults to vacuum (1/2)·I.
"""
function spectral_density(sys::QuantumStateSpace, freqs::AbstractVector{<:Real};
                           S_in = vacuum_noise(sys))
    G = freqresp(sys, freqs)           # (nout, nin, nω)
    nout, nin, nω = size(G)

    S = Array{ComplexF64}(undef, nout, nout, nω)
    for k in 1:nω
        Gk = G[:, :, k]
        S[:, :, k] = Gk * S_in * adjoint(Gk)
    end

    names = _quadrature_output_names(sys)
    return SpectralDensityMatrix(collect(Float64, freqs), S, names)
end

