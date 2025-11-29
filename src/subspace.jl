#= This file defines the Subspace type. This type's role is analogus to that of the ConcreteHilbertSpace in SecondQuantizedAlgebra.jl, and each SLH system
# has a vector of subspaces which are in one-to-one correspondence with the 'spaces' of a SecondQuantizedAlgebra.ProductSpace. Each Subspace type provides
# the names of the parameters associated with the corresponding Hilbert space. =#

abstract type Subspace end

function promote_name(mode::Subspace,parentname)
    newname = parentname*"_"*mode.name
    typeof(mode).name.wrapper(newname)
end

struct MechanicalMode <: Subspace
    name::String
end

function parameternames(subsys::MechanicalMode)
    if subsys.name == ""
        return [:Ω,:m,:Γ]
    else
        fsymb = Symbol(subsys.name,"_",:Ω)
        msymb = Symbol(subsys.name,"_",:m)
        gsymb = Symbol(subsys.name,"_",:Γ)
        return [fsymb,msymb,gsymb]
    end
end

function operatornames(subsys::MechanicalMode)
    return [:b]
end

function quadratureblocks(sys, subsys::MechanicalMode)
    params = parameternames(subsys)
    
    m_name = params[2]
    m = sys.parameters[m_name]

    w_name = params[1]
    w = sys.parameters[w_name]

    left = [0.5 0.5; -0.5im*m*w 0.5im*m*w]
    right = [1 im/(m*w); 1 -im/(m*w)]

    return (left,right)
end


struct OpticalMode <: Subspace
    name::String
end

function parameternames(subsys::OpticalMode)
    if subsys.name == ""
        return [:ω,:l,:κ]
    else
        fsymb = Symbol(subsys.name,"_",:ω)
        msymb = Symbol(subsys.name,"_",:l)
        gsymb = Symbol(subsys.name,"_",:κ)
        return [fsymb,msymb,gsymb]
    end
end

function operatornames(subsys::OpticalMode)
    return [:a]
end

function quadratureblocks(sys, subsys::OpticalMode)
    left = 1/sqrt(2)*[1 1; -im im]
    right = 1/sqrt(2)*[1 im; 1 -im]
    return (left, right) 
end

struct GenericMode <: Subspace
    name::String
end

function quadratureblocks(sys,subsys::GenericMode)
    return quadratureblocks(sys,OpticalMode(""))
end
