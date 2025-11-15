abstract type Subspace end

function promotename(mode::Subspace,parentname)
    newname = parentname*"_"*mode.name
    typeof(mode).name.wrapper(newname)
end

struct MechanicalMode <: Subspace
    name::String
end

function parameternames(subsys::MechanicalMode)
    return [:Ω,:m,:Γ]
end

function operatornames(subsys::MechanicalMode)
    return [:b]
end

struct OpticalMode <: Subspace
    name::String
end

function parameternames(subsys::OpticalMode)
    return [:ω,:L,:κ]
end

function operatornames(subsys::OpticalMode)
    return [:a]
end

struct GenericMode <: Subspace
    name::String
end
