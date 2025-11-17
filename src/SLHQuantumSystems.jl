module SLHQuantumSystems

using SecondQuantizedAlgebra
using SymbolicUtils
using LinearAlgebra
using Symbolics

include("qsymbols.jl")
export get_qnumbers, get_cnumbers, get_additive_terms, islinear, ordered_qsymbols, coeff

include("subspace.jl")
export parameternames, operatornames, quadratureblocks, OpticalMode, MechanicalMode, GenericMode

include("slh.jl")
export SLH, concatenate, feedbackreduce, operators, parameters, promote_name

import Symbolics.substitute
include("abcd.jl")
export state_vector, makedriftA, makeinputB, eqsofmotion, slh2abcd, fresponse_allIO, fresponse_state2output, toquadrature, symbfresponse, StateSpace

include("componentlibrary.jl")
export cavity, squeezing_cavity, radiation_pressure_cavity, qed_cavity

## Handle extensions
extension_fns = [
    :Makie => [:bode],
]

for (_pkg, fns) in extension_fns
    for fn in fns
        @eval function $fn end
        @eval export $fn
    end
end

function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

function __init__()
    # Hint if package extension is required to access functionality.
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, _argtypes, _kwargs
            for (pkg, fns) in extension_fns
                if in(Symbol(exc.f), fns) && !is_pkg_loaded(pkg)
                    print(io, "\nImport package ")
                    printstyled(io, "$pkg"; color = :cyan, bold = true)
                    print(io, " to enable the ")
                    printstyled(io, "$(exc.f)"; italic = true, color = :light_magenta)
                    print(io, " method.")
                end
            end
        end
    end
end


end
