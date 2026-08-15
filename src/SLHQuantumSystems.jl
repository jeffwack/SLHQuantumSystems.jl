module SLHQuantumSystems

using SecondQuantizedAlgebra
using SymbolicUtils
using LinearAlgebra
using Symbolics
import Symbolics.substitute #defined for QuantumStateSpace in abcd.jl

include("subspace.jl")
export promote_name
export parameternames, operatornames, OpticalMode, MechanicalMode, GenericMode
export quadrature_parameter_names, quadrature_transform, quadrature_scale
export zpf_length, zpf_momentum
export param_key

include("slh.jl")
export SLH, concatenate, feedbackreduce

include("abcd.jl")
export state_vector, eqsofmotion, toquadrature, QuantumStateSpace
export LadderBasis, QuadratureBasis
export ssdata, nstates, ninputs, noutputs

include("fresponse.jl")
export fresponse_state2output, symbfresponse

include("componentlibrary.jl")
export cavity, squeezing_cavity

include("spectral_density.jl")
export SpectralDensityMatrix, spectral_density, vacuum_noise

end
