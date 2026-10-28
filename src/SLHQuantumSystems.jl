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
export eqsofmotion, toquadrature, QuantumStateSpace
export LadderBasis, QuadratureBasis
export ssdata, nstates, ninputs, noutputs
export mode_names, state_names, input_ports, output_ports, input_names, output_names

include("fresponse.jl")
export TransferMatrix, transfer_matrix, symbfresponse

include("componentlibrary.jl")
export cavity, squeezing_cavity

end
