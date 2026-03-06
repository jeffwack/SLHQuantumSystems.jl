#= This file defines the Subspace type. This type's role is analogus to that of the ConcreteHilbertSpace in SecondQuantizedAlgebra.jl, and each SLH system
# has a vector of subspaces which are in one-to-one correspondence with the 'spaces' of a SecondQuantizedAlgebra.ProductSpace. Each Subspace type provides
# the names of the parameters associated with the corresponding Hilbert space. =#

using PhysicalConstants.CODATA2018: ReducedPlanckConstant as ℏ_SI

abstract type Subspace end

function promote_name(mode::Subspace,parentname)
    newname = parentname*"_"*mode.name
    typeof(mode).name.wrapper(newname)
end

"""
    param_key(mode::Subspace, base::Symbol) → Symbol

Return the parameter dict key for `base` in the context of `mode`.
For an unnamed mode (`mode.name == ""`), returns `base` directly.
For a named mode, returns `Symbol(mode.name, "_", base)`.

This is the single source of truth for parameter naming: all functions
that look up mode parameters should use this rather than positional indexing.
"""
param_key(mode::Subspace, base::Symbol) =
    mode.name == "" ? base : Symbol(mode.name, "_", base)

struct MechanicalMode <: Subspace
    name::String
end

function parameternames(subsys::MechanicalMode)
    bases = [:Ω, :m, :Γ]
    return [param_key(subsys, b) for b in bases]
end

function operatornames(subsys::MechanicalMode)
    return [:b]
end

quadrature_parameter_names(subsys::MechanicalMode) = parameternames(subsys)[1:2]

function quadrature_transform(subsys::MechanicalMode, params::Dict)
    m = params[param_key(subsys, :m)]
    w = params[param_key(subsys, :Ω)]
    left  = [0.5 0.5; -0.5im*m*w 0.5im*m*w]
    right = [1 im/(m*w); 1 -im/(m*w)]
    return (left, right)
end

"""
Zero-point fluctuation amplitude [m].
x_zpf = √(ℏ / (2 m Ω))
"""
function zpf_length(subsys::MechanicalMode, params::Dict)
    m = params[param_key(subsys, :m)]
    Ω = params[param_key(subsys, :Ω)]
    return sqrt(ℏ_SI.val / (2 * m * Ω))
end

"""Zero-point fluctuation momentum [kg⋅m/s]."""
function zpf_momentum(subsys::MechanicalMode, params::Dict)
    m = params[param_key(subsys, :m)]
    Ω = params[param_key(subsys, :Ω)]
    return sqrt(ℏ_SI.val * m * Ω / 2)
end

"""SI normalization diagonal [x_zpf, p_zpf] for the mechanical quadrature state vector."""
function quadrature_scale(subsys::MechanicalMode, params::Dict)
    return [zpf_length(subsys, params), zpf_momentum(subsys, params)]
end


struct OpticalMode <: Subspace
    name::String
end

function parameternames(subsys::OpticalMode)
    bases = [:κ, :Δ, :ω, :l]
    return [param_key(subsys, b) for b in bases]
end

function operatornames(subsys::OpticalMode)
    return [:a]
end

quadrature_parameter_names(subsys::OpticalMode) = Symbol[]

function quadrature_transform(subsys::OpticalMode, params::Dict)
    left  = 1/sqrt(2)*[1 1; -im im]
    right = 1/sqrt(2)*[1 im; 1 -im]
    return (left, right)
end

quadrature_scale(subsys::OpticalMode, params::Dict) = [1.0, 1.0]


struct GenericMode <: Subspace
    name::String
end

function parameternames(subsys::GenericMode)
    return [param_key(subsys, :ω)]
end

function operatornames(subsys::GenericMode)
    return [:a]
end

quadrature_parameter_names(subsys::GenericMode) = Symbol[]

quadrature_transform(subsys::GenericMode, params::Dict) = quadrature_transform(OpticalMode(""), params)

quadrature_scale(subsys::GenericMode, params::Dict) = [1.0, 1.0]


# Compatibility shim: keep quadratureblocks working for existing call sites
function quadratureblocks(sys, subsys::Subspace)
    relevant = Dict(k => sys.parameters[k] for k in quadrature_parameter_names(subsys))
    return quadrature_transform(subsys, relevant)
end


"""
Single-photon optomechanical coupling [rad/s].
g₀ = (ω_L / l) * x_zpf
"""
function g0_coupling(mech::MechanicalMode, opt::OpticalMode, params::Dict)
    x_zpf = zpf_length(mech, params)
    ω_L   = params[param_key(opt, :ω)]
    l     = params[param_key(opt, :l)]
    return (ω_L / l) * x_zpf
end

"""
Enhanced optomechanical coupling from intracavity photon number n̄ [rad/s].
g = g₀ * √n̄
"""
function g_coupling(mech::MechanicalMode, opt::OpticalMode, params::Dict, n_photon::Real)
    return g0_coupling(mech, opt, params) * sqrt(n_photon)
end
