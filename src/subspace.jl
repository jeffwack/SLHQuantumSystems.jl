abstract type Subspace end

function promotename(mode::Subspace,parentname)
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
    left = [0.5 0.5; -0.5im 0.5im]
    right = [1 im; 1 -im]
    return (left, right) 
end

struct GenericMode <: Subspace
    name::String
end

function quadratureblocks(sys,subsys::GenericMode)
    return quadratureblocks(sys,OpticalMode(""))
end
