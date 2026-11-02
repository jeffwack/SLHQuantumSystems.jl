# This file defines the Subspace types. Each SLH system has one subspace per factor of its Hilbert space.
# A subspace holds only a name; the names of everything found there follow from it and from two types:
#   - the type of subspace says which parameters it has, e.g. an OpticalMode has a length L
#   - the type of its Hilbert space (from SecondQuantizedAlgebra) says which operators act on it
# Each name is the subspace name and a suffix, as for states and ports: "cavity_L", "cavity_a".

"""
    Subspace{S}

A factor of an SLH system's Hilbert space, such as a mode or an atom, whose Hilbert space
has type `S` (`FockSpace`, `PhaseSpace`, `NLevelSpace`, …). Each subtype has two fields,
`name::String` and `spaceargs::Tuple`, the arguments of `S` after its name, such as the
number of levels of an `NLevelSpace`.
"""
abstract type Subspace{S} end

"""
    GenericSubspace(name, space, spaceargs...)

A subspace with no parameters, on a Hilbert space of type `space`:
`GenericSubspace("atom", NLevelSpace, 2)`.
"""
struct GenericSubspace{S} <: Subspace{S}
    name::String
    spaceargs::Tuple
end
# These constructors take a name from the user and check it; promote_name builds promoted names directly
GenericSubspace(name, space::Type, spaceargs...) = (check_name("subspace", name); GenericSubspace{space}(name, spaceargs))

"""
    OpticalMode(name, space = FockSpace, spaceargs...)

An optical mode, such as the mode of a cavity. Its parameter is the length `L` [m].
"""
struct OpticalMode{S} <: Subspace{S}
    name::String
    spaceargs::Tuple
end
OpticalMode(name, space::Type = FockSpace, spaceargs...) = (check_name("subspace", name); OpticalMode{space}(name, spaceargs))

"""
    MechanicalMode(name, space = FockSpace, spaceargs...)

A mechanical mode, such as the motion of a mirror. Its parameter is the mass `m` [kg].
`MechanicalMode("mirror", PhaseSpace)` has position and momentum operators, as for a free mass.
"""
struct MechanicalMode{S} <: Subspace{S}
    name::String
    spaceargs::Tuple
end
MechanicalMode(name, space::Type = FockSpace, spaceargs...) = (check_name("subspace", name); MechanicalMode{space}(name, spaceargs))

# The subspace's Hilbert space, built here rather than read from a ProductSpace
hilbertspace(sub::Subspace{S}) where {S} = S(Symbol(sub.name), sub.spaceargs...)

# The name is promoted as in S, L and H, and the names found at the subspace follow it
function promote_name(sub::Subspace, parentname)
    newname = isempty(sub.name) ? parentname : parentname * "_" * sub.name
    return typeof(sub)(newname, sub.spaceargs)
end

# The name of something found at a subspace: "cavity_L", or "L" for an unnamed subspace
subname(sub, suffix) = Symbol(isempty(sub.name) ? suffix : sub.name * "_" * suffix)

# The parameters found at a subspace, by its type
parameter_suffixes(::GenericSubspace) = ()
parameter_suffixes(::OpticalMode) = ("L",)
parameter_suffixes(::MechanicalMode) = ("m",)
subspace_parameters(sub) = [Symbolics.variable(subname(sub, s)) for s in parameter_suffixes(sub)]

# The operators found at a subspace, by the type of its Hilbert space. `h` is the system's Hilbert
# space and `at` the subspace's index in it, empty when the system has only this subspace.
subspace_operators(sub::Subspace{FockSpace}, h, at...) = [Destroy(h, subname(sub, "a"), at...)]
function subspace_operators(sub::Subspace{PhaseSpace}, h, at...)
    return [Position(h, subname(sub, "x"), at...), Momentum(h, subname(sub, "p"), at...)]
end
function subspace_operators(sub::Subspace{NLevelSpace}, h, at...)
    levels = first(sub.spaceargs)
    levels = levels isa Integer ? (1:levels) : levels
    return [Transition(h, subname(sub, "σ"), i, j, at...) for i in levels, j in levels]
end

#rename keeps the symtype and metadata, so a complex (::Number) parameter stays complex
function promote_name(parameter, parentname)
    old_sym = nameof(parameter)
    new_sym = Symbol(parentname,"_",old_sym)
    return Symbolics.rename(parameter, new_sym)
end
